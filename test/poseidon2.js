import poseidon2Constants from "./poseidon2_constants.js";
import { utils } from "ffjavascript";
import { Poseidon2, F1Field } from "poseidon2";

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
//ex: For Vesta prime, t=3 https://github.com/HorizenLabs/poseidon2/blob/main/plain_implementations/src/poseidon2/poseidon2_instance_vesta.rs

const field = new F1Field(
    BigInt(
        "21888242871839275222246405745257275088548364400416034343698204186575808495617"
    )
);

const instance = new Poseidon2(
    getPoseidon2Params(3, 5, 8, 56, MAT_DIAG3_M_1, MAT_INTERNAL3, RC3),
    field
);

const state = [BigInt(0), BigInt(1), BigInt(2)];
const permuteResult = instance.permute(state);
console.log(permuteResult);