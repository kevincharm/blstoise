import { expect } from 'chai'
import { pair, validatePairing } from '../src/pairing'
import { PointG1, PointG2 } from '../src/point'
import { Fq, Fq12, Fq2 } from '../src/ff'
import testVectors from './vectors/bls_pairing'
import { computeWitness, verifyEquivalentPairings } from '../src/witness'

function toFp(hex: Buffer): Fq {
    return new Fq(BigInt(`0x${hex.toString('hex')}`))
}

describe('pairing', () => {
    const g1 = PointG1.generator()
    const g2 = PointG2.generator()

    it('neg G1', () => {
        const x = pair(g1, g2)
        const y = pair(g1.neg(), g2)
        expect(x.mul(y).equals(Fq12.one())).to.eq(true)
        expect(y.equals(x)).to.eq(false, 'pairing is non-degenerate')
    })

    it('neg G2', () => {
        const x = pair(g1, g2)
        const y = pair(g1, g2.neg())
        expect(x.mul(y).equals(Fq12.one())).to.eq(true)
        expect(y.equals(x)).to.eq(false, 'pairing is non-degenerate')
    })

    it('bilinearity in G1', () => {
        const x = pair(g1, g2)
        const y = pair(g1.mul(2n), g2)
        expect(x.mul(x).equals(y)).to.eq(true)
        expect(y.equals(x)).to.eq(false, 'pairing is non-degenerate')
    })

    it('bilinearity in G2', () => {
        const x = pair(g1, g2)
        const y = pair(g1, g2.mul(2n))
        expect(x.mul(x).equals(y)).to.eq(true)
        expect(y.equals(x)).to.eq(false, 'pairing is non-degenerate')
    })

    it('composite', () => {
        const x = pair(g1.mul(25n), g2.mul(42n))
        const y = pair(g1.mul(1050n), g2)
        expect(x.equals(y)).to.eq(true)
    })

    for (const testVector of testVectors) {
        const { ps, qs, result } = parseGethTest(testVector)

        it(`[geth] ${testVector.Name} (${ps.length} pairings)`, () => {
            expect(validatePairing(ps, qs)).to.eq(result)
        }).timeout(0) // pairing is mad slow bruv

        // Test computing witness residues for valid pairings
        if (result) {
            it(`[geth] ${testVector.Name} (${ps.length} pairings) -> compute witness residues`, () => {
                let f = Fq12.one()
                for (let i = 0; i < ps.length; i++) {
                    f = f.mul(pair(ps[i], qs[i]))
                }
                const { c, wi } = computeWitness(f)
                expect(verifyEquivalentPairings(ps, qs, c, wi)).to.eq(true)
            }).timeout(0)
        }
    }

    it('compute and verify witness residues', () => {
        // Get equivalent pairings x and y
        const ps = [g1, g1.neg()]
        const qs = [g2, g2]
        const x = pair(g1, g2)
        const y = pair(g1.neg(), g2)
        // Compute the witness for equivalent pairing check
        const { c, wi } = computeWitness(x.mul(y))
        expect(verifyEquivalentPairings(ps, qs, c, wi)).to.eq(true)
    }).timeout(0)
})

function parseGethTest(vector: (typeof testVectors)[number]) {
    const input = Buffer.from(vector.Input, 'hex')
    if (input.byteLength % 384 !== 0) {
        throw new Error(`Invalid input length: ${vector.Name}`)
    }
    const len = input.byteLength / 384
    const result = Boolean(BigInt(`0x${vector.Expected}`))
    const ps: PointG1[] = []
    const qs: PointG2[] = []
    for (let i = 0; i < len; i++) {
        const cur = i * 384
        const p = new PointG1(
            toFp(input.subarray(cur, cur + 64)),
            toFp(input.subarray(cur + 64, cur + 128)),
        )
        const q = new PointG2(
            new Fq2(
                toFp(input.subarray(cur + 128, cur + 192)),
                toFp(input.subarray(cur + 192, cur + 256)),
            ),
            new Fq2(
                toFp(input.subarray(cur + 256, cur + 320)),
                toFp(input.subarray(cur + 320, cur + 384)),
            ),
        )
        ps.push(p)
        qs.push(q)
    }
    return {
        ps,
        qs,
        result,
    }
}
