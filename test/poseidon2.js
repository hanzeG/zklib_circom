const chai = require("chai");
const path = require("path");
const wasm_tester = require("circom_tester").wasm;
const assert = chai.assert;

const poseidon2Constants = require("./poseidon2_constants.js");
const { utils, getCurveFromName } = require("ffjavascript");
const { Poseidon2, F1Field } = require("poseidon2");

const { MAT_DIAG3_M_1, MAT_INTERNAL3, RC3 } = utils.unstringifyBigInts(poseidon2Constants);

// call this function with your parameters from sage/horizen labs' precomputed constants
function getPoseidon2Params(
    t,
    d,
    rounds_f,
    rounds_p,
    mat_internal_diag_m_1,
    mat_internal,
    round_constants
) {
    const r = rounds_f / 2;
    const rounds = rounds_f + rounds_p;
    return {
        t: t,
        d: d,
        rounds_f_beginning: r,
        rounds_p: rounds_p,
        rounds_f_end: r,
        rounds: rounds,
        mat_internal_diag_m_1: mat_internal_diag_m_1,
        _mat_internal: mat_internal,
        round_constants: round_constants,
    };
}

describe("Poseidon2 Circuit test", function () {
    let prime;
    let F;
    let poseidon2;
    let circuit;

    this.timeout(1000000);

    before(async () => {
        prime = await getCurveFromName("bn128", true);
        F = new F1Field(prime.r);
        poseidon2 = new Poseidon2(
            getPoseidon2Params(3, 5, 8, 56, MAT_DIAG3_M_1, MAT_INTERNAL3, RC3), F
        );
        circuit = await wasm_tester(path.join(__dirname, "circuits", "poseidon2_3_test.circom"));
    });

    it("Should check constrain of hash([1, 1, 1]) t=3", async () => {
        const input = [1, 1, 1];
        const bigIntInput = input.map(x => BigInt(x));
        const res2 = poseidon2.permute(bigIntInput);

        const w = await circuit.calculateWitness({ inputs: input });

        await circuit.assertOut(w, { out: res2 });
        await circuit.checkConstraints(w);
    });

});