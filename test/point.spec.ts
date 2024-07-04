import { expect } from 'chai'
import { PointG1, PointG2 } from '../src/point'
import { randomBigInt } from '../src/utils'
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

        it.skip('add', () => {
            const p = PointG1.generator()
            const p2 = p.add(p)
            const p3 = p2.add(p)
            // TODO: Test against ref
        })

        it('mul', () => {
            const p = PointG1.generator()
            const p2 = p.mul(2n)
            const p2Add = p.add(p)
            expect(p2.equals(p2Add)).to.eq(true)
            const p3 = p.mul(3n)
            expect(p3.equals(p2Add.add(p))).to.eq(true)
        })

        it('is on curve / subgroup check', () => {
            const p = PointG1.generator().mul(randomBigInt(31))
            expect(p.isOnCurve()).to.eq(true)
            expect(p.isInSubgroup()).to.eq(true)

            const pWrong = new PointG1(new Fq(1n), new Fq(0n))
            expect(pWrong.isOnCurve()).to.eq(false)
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

        it('is on curve / subgroup check', () => {
            const p = PointG2.generator().mul(randomBigInt(31))
            expect(p.isOnCurve()).to.eq(true)
            expect(p.isInSubgroup()).to.eq(true)

            const pWrong = new PointG2(Fq2.fromTuple([1n, 0n]), Fq2.fromTuple([0n, 0n]))
            expect(pWrong.isOnCurve()).to.eq(false)
        })
    })
})
