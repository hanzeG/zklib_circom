pragma circom 2.0.0;

include "./poseidon2_constants.circom";

template Sigma() {
    signal input in;
    signal output out;

    signal in2;
    signal in4;

    in2 <== in*in;
    in4 <== in2*in2;

    out <== in4*in;
}

template Rc(t, RC, r) {
    signal input in[t];
    signal output out[t];

    for (var i=0; i<t; i++) {
        out[i] <== in[i] + RC[r][i];
    }
}

template Sum(t) {
    signal input in[t];
    signal output out;
    var sum = 0;
    for (var i = 0; i < t; i++) {
        sum += in[i];
    }
    out <== sum;
}

template Internal(t) {
    signal input in[t];
    signal output out[t];

    var diagMat[t] = POSEIDON2_DIAG(t);

    component sum = Sum(t);

    for (var i = 0; i < t; i++) {
        sum.in[i] <== in[i];
    }
    for (var i = 0; i < t; i++) {
        out[i] <== sum.out + in[i] * diagMat[i];
    }
}

template Matmul_m4(t) {
    signal input in[t];
    signal output out[t];
    
    // [[5, 7, 1, 3]
    // [4, 6, 1, 1]
    // [1, 3, 5, 7]
    // [1, 1, 4, 6]]
    var t4 = t/4;
    for (var i =0; i < t4; i++){        
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
    signal input in[t];
    signal output out[t];
    
    component sum;
    component m4;

    if (t < 4) {
        sum = Sum(t);
        for (var i = 0; i < t; i++) {
            sum.in[i] <== in[i];
        }
        for (var i = 0; i < t; i++) {
            out[i] <== in[i] + sum.out;
        }
    } else if (t == 4) {
        m4 = Matmul_m4(t);
        for (var i = 0; i < t; i++) {
            m4.in[i] <== in[i];
        }
        for (var i = 0; i < t; i++) {
            out[i] <== m4.out[i];
        }
    } else { // t > 4
        m4 = Matmul_m4(t);
        var tmp[t];
        for (var i = 0; i < t; i++) {
            m4.in[i] <== in[i];
        }
        for (var i = 0; i < t; i++) {
            tmp[i] = m4.out[i];
        }

        var t4 = t / 4;
        var store[4];
        for (var i = 0; i < 4; i++) {
            store[i] = tmp[i];
            for (var j = 1; j < t4; j++) {
                store[i] += tmp[i + 4 * j];
            }
        }

        for (var i = 0; i < t4; i++) {
            out[i * 4 + 0] <== store[0] + in[i * 4 + 0];
            out[i * 4 + 1] <== store[1] + in[i * 4 + 1];
            out[i * 4 + 2] <== store[2] + in[i * 4 + 2];
            out[i * 4 + 3] <== store[3] + in[i * 4 + 3];
        }
    }
}


template Poseidon2Ex(nInputs) {
    signal input inputs[nInputs];
    signal output out[nInputs];

    var t = nInputs;
    var nRoundsF = 8;
    var N_ROUNDS_P = 56;
    var nRoundsP = N_ROUNDS_P;
    var RC[nRoundsP+nRoundsF][t] = POSEIDON2_RC(t);

    component rc[nRoundsF+nRoundsP];
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
    for (var i=0; i< nRoundsF/2; i++) {
        rc[i] = Rc(t, RC, i);
        external[i] = External(t);
        for (var j=0; j<t; j++) {
            rc[i].in[j] <== currentState[j];
        }
        for (var j=0; j<t; j++) {
            sigmaF[i][j] = Sigma();
            sigmaF[i][j].in <== rc[i].out[j];
            external[i].in[j] <== sigmaF[i][j].out;
        }
        for (var j=0; j<t; j++) {
            currentState[j] = external[i].out[j];
        }
    }

    // internal round 2
    for (var i=nRoundsF/2; i<nRoundsF/2+nRoundsP; i++) {
        rc[i] = Rc(1, RC, i);
        sigmaP[i-nRoundsF/2] = Sigma();
        internal[i-nRoundsF/2] = Internal(t);

        rc[i].in[0] <== currentState[0];
        sigmaP[i-nRoundsF/2].in <== rc[i].out[0];
        currentState[0] = sigmaP[i-nRoundsF/2].out;
        
        for (var j=0; j<t; j++) {
            internal[i-nRoundsF/2].in[j] <== currentState[j];
        }
        for (var j=0; j<t; j++) {
            currentState[j] = internal[i-nRoundsF/2].out[j];
        }
    }

    // exteranl round 3
    for (var i=nRoundsF/2+nRoundsP; i<nRoundsF+nRoundsP; i++) {
        rc[i] = Rc(t, RC, i);
        external[i-nRoundsP] = External(t);
        for (var j=0; j<t; j++) {
            rc[i].in[j] <== currentState[j];
        }
        for (var j=0; j<t; j++) {
            sigmaF[i-nRoundsP][j] = Sigma();
            sigmaF[i-nRoundsP][j].in <== rc[i].out[j];
            external[i-nRoundsP].in[j] <== sigmaF[i-nRoundsP][j].out;
        }
        for (var j=0; j<t; j++) {
            currentState[j] = external[i-nRoundsP].out[j];
        }
    }

    for (var i=0; i<nInputs; i++) {
        out[i] <== currentState[i];
    }
}

template Poseidon2(nInputs,nOutputs) {
    signal input inputs[nInputs];
    signal output out[nOutputs];

    component pEx = Poseidon2Ex(nInputs);
    for (var i=0; i<nInputs; i++) {
        pEx.inputs[i] <== inputs[i];
    }
    for (var i=0; i<nOutputs; i++) {
        out[i] <== pEx.out[i];
    }
}
