import { abs, Fq12, X } from './ff'
import { PointG1, PointG2 } from './point'

// Compute asymmetric pairing e(p,q)
export function pair(p: PointG1, q: PointG2): Fq12 {
    // Base case: p or q is the point at infinity
    if (p.equals(PointG1.infinity()) || q.equals(PointG2.infinity())) {
        return Fq12.one()
    }
    if (!p.isInSubgroup() || !p.isOnCurve()) {
        throw new Error(`Invalid point: ${p}`)
    }
    if (!q.isInSubgroup() || !q.isOnCurve()) {
        throw new Error(`Invalid point: ${q}`)
    }
    const r = miller(p, q)
    return r.finalExp()
}

export function validatePairing(ps: PointG1[], qs: PointG2[]): boolean {
    if (ps.length !== qs.length) {
        throw new Error('Invalid inputs')
    }

    let result = Fq12.one()
    for (let i = 0; i < ps.length; i++) {
        const p = ps[i]
        const q = qs[i]
        if (!p.isOnCurve() || !p.isInSubgroup()) throw new Error(`Invalid point: ${p}`)
        if (!q.isOnCurve() || !q.isInSubgroup()) throw new Error(`Invalid point: ${q}`)
        const r = miller(p, q)
        result = result.mul(r)
    }
    result = result.finalExp()
    return result.equals(Fq12.one())
}

export function miller(p: PointG1, q: PointG2): Fq12 {
    // Binary representation of curve parameter B
    // NB: This can be precomputed!
    const iterations: boolean[] = []
    let curveX = abs(X)
    while (curveX > 0n) {
        const isOddBit = Boolean(curveX & 1n)
        iterations.push(isOddBit)
        curveX >>= 1n
    }
    iterations.pop()
    iterations.reverse()

    // Miller loop
    let acc = Fq12.one()
    let r = q.clone()
    for (const i of iterations) {
        const doubleR = r.double()
        acc = acc.mul(acc).mul(lineDouble(r, p))

        if (i) {
            r = doubleR.add(q)
            acc = acc.mul(lineAdd(doubleR, q, p))
        } else {
            r = doubleR
        }
    }
    return acc
}

function lineDouble(r: PointG2, p: PointG1): Fq12 {
    const wideR = r.untwist()
    const slope = wideR.x
        .mul(wideR.x)
        .mul(Fq12.fromNumber(3n))
        .mul(wideR.y.mul(Fq12.fromNumber(2n)).inv())
    const v = wideR.y.sub(slope.mul(wideR.x))
    return Fq12.fromNumber(p.y.value).sub(Fq12.fromNumber(p.x.value).mul(slope)).sub(v)
}

function lineAdd(r: PointG2, q: PointG2, p: PointG1): Fq12 {
    const wideR = r.untwist()
    const wideQ = q.untwist()
    if (wideR.x.equals(wideQ.x) && wideR.y.equals(wideQ.y.neg())) {
        return Fq12.fromNumber(p.x.value).sub(wideR.x)
    } else {
        const slope = wideQ.y.sub(wideR.y).mul(wideQ.x.sub(wideR.x).inv())
        const v = wideQ.y.mul(wideR.x).sub(wideR.y.mul(wideQ.x)).mul(wideR.x.sub(wideQ.x).inv())
        return Fq12.fromNumber(p.y.value).sub(Fq12.fromNumber(p.x.value).mul(slope)).sub(v)
    }
}
