pragma circom 2.1.9;

include "./poseidon2_constants.circom";
include "../node_modules/circomlib/circuits/poseidon.circom"

template SumMat(t) {
    signal input in[t];
    signal output out;
    var sum = 0;
    for (var i = 0; i < t; i++) {
        sum <== sum + in[i];
    }
    out <== sum;
}

template DotMat(t, MDS) {
    signal input in[t];
    signal output out[t];

    for (var i = 0; i < t; i++) {
        out[i] <== 0;
        for (var j = 0; j < n; j++) {
            out[i] <== out[i] + MDS[i][j] * in[j];
        }
    }
}

template DiagMat(t) {
    signal input in[t][t];
    signal output out[t];

    for (var i = 0; i < t; i++) {
        out[i] <== in[i][i];
    }
}

template Internal(t, MDS) {
    signal input in[t];
    signal output out[t];

    if (t < 4) {
        component dot = Dot(t);
        for (var i = 0, i < t, i++) {
            dot.inA[i] <== in[i];
            dot.inB[i] <== MDS[i];
        }
        for (var i = 0, i < t, i++) {
            out[i] <== dot.out[i];
        }
    } else {
        var diagMat[t];
        component diag = Diag(t);
        component sum = Sum(t);
        for (var i = 0, i < t, i++) {
            for (var j = 0, j < t, j++) {
                diag.in[i][j] <== MDS[i][j];
            }
            sum.in[i] <== in[i];
        }
        for (var i = 0; i < t; i++) {
            out[i] <== in[i] * diag.out[i] - in[i] + sum.out;
        }
    }
}

template Matmul_m4(t) {
    input signal in[t];
    output signal out[t];
    
    var t4 = t/4;
    for (var i =0, i < t4, i++){        
        var tmp[8];
        tmp[0] = in[0+i*4] + in[1+i*4];
        tmp[1] = in[2+i*4] + in[3+i*4];
        tmp[2] = in[1+i*4] * 2 + tmp[1];
        tmp[3] = in[3+i*4] * 2 + tmp[0];
        tmp[4] = tmp[1] * 4 + tmp[3];
        tmp[5] = tmp[0] * 4 + tmp[2];
        tmp[6] = tmp[3] + tmp[5];
        tmp[7] = tmp[2] + tmp[4];

        out[0+i*4] <== tmp[6];
        out[1+i*4] <== tmp[5];
        out[2+i*4] <== tmp[7];
        out[3+i*4] <== tmp[4];
    }
}

template External(t) {
    input signal in[t];
    input signal out[t];

    var tmp[t];

    if (t < 4) {
        component sum = SumMat(t);
        for (var i = 0, i < t, i++) {
            sum.in[i] <== in[i];
        }
        for (var i = 0; i < t; i++) {
            out[i] <== in[i] + sum;
        }
    } else if (t >= 4) {
        component m4 = Matmul_m4(t);

        for (var i = 0, i < t, i++) {
            m4.in[i] <== in[i];
        }
        for (var i = 0, i < t, i++) {
            tmp[i] <== m4.out[i];
        }

        if (t > 4) {
            var t4 = t/4;
            var store[4];
            for (var i=0, i<4, i++){
                store[i] = tmp[i];
                for (var j=1, j<t4, j++){
                store[i] += tmp[i+4*j];
                }
            }

            for (var i=0, i<t4, i++){
                out[i*4+0] <== store[0] + in[i*4+0];
                out[i*4+1] <== store[1] + in[i*4+1];
                out[i*4+2] <== store[2] + in[i*4+2];
                out[i*4+3] <== store[3] + in[i*4+3];
            }
        } else if {
            for (var i = 0, i < t, i++) {
                out[i] <== tmp[i];
            }
        }        
    }
}

template Poseidon2Ex(nInputs, nOuts) {
    signal input inputs[nInputs];
    signal output out[nOuts];

    var t = nInputs;
    var nRoundsF = 8;
    var N_ROUNDS_P = 56;
    var nRoundsP = N_ROUNDS_P;
    var C_INT[nRoundsP] = POSEIDON2_C_INT(t);
    var C_EXT[t][nRoundsF] = POSEIDON2_C_EXT(t);

    component arkF[nRoundsF];
    component arkP[nRoundsP];
    component sigmaF[nRoundsF][t];
    component sigmaP[nRoundsP];
    component initialEx;
    component external[nRoundsF];
    component internal[nRoundsP];

    // Linear layer at beginning
    initialEx = External(t);
    for (var i=0; i<t; i++) {
        initialEx.in[i] <== inputs[i];
    }
    var currentState[t];
    for (var i=0; i<t; i++) {
        currentState[i] = initialEx.out[i];
    }

    // exteranl round 1
    for (var i=0; i<nRoundsF/2; i++) {
        arkF[i] = Ark[t, C_INT, i*t];
        external[i] = External(t);
        for (var j=0; j<t; j++) {
            arkF[i].in[j] <== currentState[j];
        }
        for (var j=0; j<t; j++) {
            sigmaF[i][j] = Sigma();
            sigmaF[i][j].in <== arkF[i].out[j];
            external[j] <== sigmaF[i][j].out;
        }
        for (var j=0; j<t; j++) {
            currentState[j] <== external[i].out[j];
        }
    }

    // internal round 2
    for (var i=0; i<nRoundsP; i++) {
        arkP[i] = Ark(1, C_INT, (nRoundsF/2+i)*t);
        sigmaP[i] = Sigma;
        internal[i] = Internal(t, C_EXT);

        arkP[i].in[0] <== currentState[0];
        sigmaP[i].in <== arkP[i].out[0];
        currentState[0] = sigmaP[i].out;
        
        for (var j=0; j<t; j++) {
            internal[i].in[j] <== currentState[j];
        }
        for (var j=0; j<t; j++) {
            currentState[j] = internal[i].out[j];
        }
    }

    // exteranl round 3
    for (var i=0; i<nRoundsF/2; i++) {
        arkF[nRoundsF/2+i] = Ark[t, C_INT, (nRoundsF/2+nRoundsP+i)*t];
        external[nRoundsF/2+i] = External(t);
        for (var j=0; j<t; j++) {
            arkF[nRoundsF/2+i].in[j] <== currentState[j];
        }
        for (var j=0; j<t; j++) {
            sigmaF[nRoundsF/2+i][j] = Sigma();
            sigmaF[nRoundsF/2+i][j].in <== arkF[nRoundsF/2+i].out[j];
            external[j] <== sigmaF[nRoundsF/2+i][j].out;
        }
        for (var j=0; j<t; j++) {
            currentState[j] <== external[i].out[j];
        }
        for (var j=0; j<nOuts; j++) {
            out[j] <== currentState[j];
        }
    }
}

template Poseidon2(nInputs) {
    signal input inputs[nInputs];
    signal output out;

    component pEx = Poseidon2Ex(nInputs, 1);
    for (var i=0; i<nInputs; i++) {
        pEx.inputs[i] <== inputs[i];
    }
    out <== pEx.out[0];
}
