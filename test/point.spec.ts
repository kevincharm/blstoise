import { expect } from 'chai'
import { PointG1, PointG2 } from '../src/point'
import { randomBigInt, randomFq, randomFq2 } from '../src/utils'
import { Fq, Fq2 } from '../src/ff'

describe('subgroup points', () => {
    describe('G1', () => {
        it('point at infinity', () => {
            const p = PointG1.infinity()
            expect(p.add(p)).to.eq(p) // O + O = O
            const pg = PointG1.generator()
            expect(p.add(pg)).to.eq(pg) // O + G = G
            expect(pg.add(p)).to.eq(pg) // G + O = G
        })

        it('add', () => {
            const p = PointG1.generator()
            const p2 = p.add(p).add(p)
            expect(p2.equals(p.mul(3n))).to.eq(true)
        })

        it('mul', () => {
            const p = PointG1.generator()
            const p2 = p.mul(2n)
            const p2Add = p.add(p)
            expect(p2.equals(p2Add)).to.eq(true)
            const p3 = p.mul(3n)
            expect(p3.equals(p2Add.add(p))).to.eq(true)
        })

        it('is on curve', () => {
            const p = PointG1.generator().mul(randomBigInt(31))
            expect(p.isOnCurve()).to.eq(true)
            expect(p.isInSubgroup()).to.eq(true)

            const pWrong = new PointG1(new Fq(1n), new Fq(0n))
            expect(pWrong.isOnCurve()).to.eq(false)
        })

        it('subgroup check', () => {
            // Generate x until we find a valid y
            let x: Fq
            let y: Fq
            for (;;) {
                // y^2 = x^3 + 4
                x = randomFq()
                const rhs = x.mul(x).mul(x).add(new Fq(4n))
                try {
                    y = rhs.sqrt()
                    break
                } catch (_) {
                    continue
                }
            }
            const p = new PointG1(x, y)
            expect(p.isOnCurve()).to.eq(true)
            expect(p.isInSubgroup()).to.eq(false)

            const pValid = p.clearCofactor()
            expect(pValid.isOnCurve()).to.eq(true)
            expect(pValid.isInSubgroup()).to.eq(true)
        })
    })

    describe('G2', () => {
        it('point at infinity', () => {
            const p = PointG2.infinity()
            expect(p.add(p)).to.eq(p) // O + O = O
            const pg = PointG2.generator()
            expect(p.add(pg)).to.eq(pg) // O + G = G
            expect(pg.add(p)).to.eq(pg) // G + O = G
        })

        it('add', () => {
            const p = PointG2.generator()
            const p2 = p.add(p)
            const p3 = p2.add(p)
            // TODO: Test against ref
        })

        it('mul', () => {
            const p = PointG2.generator()
            const p2 = p.mul(2n)
            const p2Add = p.add(p)
            expect(p2.equals(p2Add)).to.eq(true)
            const p3 = p.mul(3n)
            expect(p3.equals(p2Add.add(p))).to.eq(true)
        })

        it('is on curve', () => {
            const p = PointG2.generator().mul(randomBigInt(31))
            expect(p.isOnCurve()).to.eq(true)
            expect(p.isInSubgroup()).to.eq(true)

            const pWrong = new PointG2(Fq2.fromTuple([1n, 0n]), Fq2.fromTuple([0n, 0n]))
            expect(pWrong.isOnCurve()).to.eq(false)
        })

        it('subgroup check', () => {
            // Generate x until we find a valid y
            let x: Fq2
            let y: Fq2
            for (;;) {
                // y^2 = x^3 + 4
                x = randomFq2()
                const rhs = x
                    .mul(x)
                    .mul(x)
                    .add(Fq2.fromTuple([4n, 4n]))
                try {
                    y = rhs.sqrt()
                    break
                } catch (_) {
                    continue
                }
            }
            const p = new PointG2(x, y)
            expect(p.isOnCurve()).to.eq(true)
            expect(p.isInSubgroup()).to.eq(false)

            const pValid = p.clearCofactor()
            expect(pValid.isOnCurve()).to.eq(true)
            expect(pValid.isInSubgroup()).to.eq(true)
        })
    })
})
