// From Novakovic & Eagen, "On Proving Pairings" https://eprint.iacr.org/2024/640.pdf
// Ref bls-lambda-roots.ipynb from Novakovic: https://github.com/akinovak/garaga/blob/4af45d5fd88bc92bfb775cd1c5c8836edfcd5c68/hydra/hints/bls-lambda-roots.ipynb
// Ref tower_to_direct_extension.py from feltroidprime: https://gist.github.com/feltroidprime/bd31ab8e0cbc0bf8cd952c8b8ed55bf5
//  and "Faster Extension Field multiplications for Emulated Pairing Circuits": https://hackmd.io/@feltroidprime/B1eyHHXNT
import { egcd, Fq12, mod } from './ff'
import { assert } from './utils'

const x = -0xd201000000010000n
const k = 12n
const r = x ** 4n - x ** 2n + 1n
const q = ((x - 1n) ** 2n / 3n) * r + x
const h = (q ** k - 1n) / r

const lambda =
    4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129030796414117214202539n
const m = lambda / r

const p = 5044125407647214251n
const h3 =
    2366356426548243601069753987687709088104621721678962410379583120840019275952471579477684846670499039076873213559162845121989217658133790336552276567078487633052653005423051750848782286407340332979263075575489766963251914185767058009683318020965829271737924625612375201545022326908440428522712877494557944965298566001441468676802477524234094954960009227631543471415676620753242466901942121887152806837594306028649150255258504417829961387165043999299071444887652375514277477719817175923289019181393803729926249507024121957184340179467502106891835144220611408665090353102353194448552304429530104218473070114105759487413726485729058069746063140422361472585604626055492939586602274983146215294625774144156395553405525711143696689756441298365274341189385646499074862712688473936093315628166094221735056483459332831845007196600723053356837526749543765815988577005929923802636375670820616189737737304893769679803809426304143627363860243558537831172903494450556755190448279875942974830469855835666815454271389438587399739607656399812689280234103023464545891697941661992848552456326290792224091557256350095392859243101357349751064730561345062266850238821755009430903520645523345000326783803935359711318798844368754833295302563158150573540616830138810935344206231367357992991289265295323280n

assert(h === 27n * p * h3)
assert(m === 3n * p ** 2n)

assert(egcd(3n, h3)[0] === 1n)
assert(egcd(p ** 2n, h3)[0] === 1n)
assert(egcd(3n, h3)[0] === 1n)
assert(egcd(p, h3)[0] === 1n)
assert(egcd(p, 27n * h3)[0] === 1n)
assert(egcd(27n, p * h3)[0] === 1n)

const ONE = Fq12.one()

/** 27th-root of unity, from sage:
 *
 * ```sage
 * q = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
 * F12.<w> = GF(q^12, modulus=x^12-2*x^6+2)
 * w27 = F12(1).nth_root(27)
 *
 * # Fq12 direct repr (sage, py_ecc) to tower repr bijection for BLS12-381
 * def direct_to_tower(x: list):
 *    res = 12 * [0]
 *    res[0] = (x[0] + x[6]) % q
 *    res[1] = x[6]
 *    res[2] = (x[2] + x[8]) % q
 *    res[3] = x[8]
 *    res[4] = (x[4] + x[10]) % q
 *    res[5] = x[10]
 *    res[6] = (x[1] + x[7]) % q
 *    res[7] = x[7]
 *    res[8] = (x[3] + x[9]) % q
 *    res[9] = x[9]
 *    res[10] = (x[5] + x[11]) % q
 *    res[11] = x[11]
 *    return res
 *
 * direct_to_tower(w27)
 * ```
 */
const w27 = Fq12.fromTuple([
    [
        [0n, 0n],
        [0n, 0n],
        [
            0n,
            3579939510681526869890171347332239103550277585524001119547067858135937309331998144326574624416802638795006015528254n,
        ],
    ],
    [
        [0n, 0n],
        [0n, 0n],
        [0n, 0n],
    ],
])
assert(w27.exp(27n).equals(ONE))
assert(!w27.exp(9n).equals(ONE))

function modinv(a: bigint, m: bigint): bigint {
    const [g, x, y] = egcd(a, m)
    if (g !== 1n) {
        throw new Error(`${a}^{-1} (mod ${m}) d.n.e.`)
    } else {
        return mod(x, m)
    }
}

function isPthResidue(x: Fq12) {
    return x.exp(h3 * 27n).equals(ONE)
}

function getPthRootInverse(x: Fq12) {
    if (isPthResidue(x)) {
        return ONE
    }

    const v = 27n * h3
    const wj = x.exp(v)

    const v_inv = modinv(v, p)
    const s = mod(-v_inv, p)

    return wj.exp(s)
}

function getOrderOf3rdPrimitiveRoot(x: Fq12) {
    // correct way is do do r * p * h3 but outputs of equal Tate pairings are always of the form c^r thus there is no rth root contribution
    const y = x.exp(p * h3)

    if (y.equals(ONE)) return 0n

    if (y.exp(3n).equals(ONE)) return 1n

    if (y.exp(9n).equals(ONE)) return 2n

    if (y.exp(27n).equals(ONE)) return 3n

    throw new Error(`Invalid ${x}`)
}

function getAny27thRootInverse(x: Fq12) {
    const pw = getOrderOf3rdPrimitiveRoot(x)

    if (pw === 0n) {
        return ONE
    }

    const _ord = 3n ** pw

    const v = p * h3
    const wj = x.exp(v)

    const v_inv = modinv(v, _ord)
    const s = mod(-v_inv, _ord)

    return wj.exp(s)
}

function h3OrdElementLambdaRoot(x: Fq12) {
    // after appying shifts we know that element is order just h3
    const e = modinv(lambda, h3)
    return x.exp(e)
}

export function computeWitness(f: Fq12) {
    // `f` will always be an r-th residue in the honest execution
    const wp_shift = getPthRootInverse(f)
    const w27_shift = getAny27thRootInverse(f)
    // witness scaling w_i
    const wi = wp_shift.mul(w27_shift)
    const f_shifted = f.mul(wi)

    const c = h3OrdElementLambdaRoot(f_shifted)
    assert(f_shifted.equals(c.exp(lambda)), 'Computation failed')

    // Actual witness residue check
    // We have (c, w^i) s.t. c^λ = f * w^i
    const c_inv = c.inv()
    assert(c_inv.exp(lambda).mul(f).mul(wi).equals(Fq12.one()), 'Witness residue check failed')

    return {
        c,
        wi,
    }
}

export function verifyEquivalentPairings(fs: [Fq12, Fq12], c: Fq12, wi: Fq12): boolean {
    const [p, q] = fs
    const f = p.mul(q)
    const c_inv = c.inv()
    return c_inv.exp(lambda).mul(f).mul(wi).equals(Fq12.one())
}
