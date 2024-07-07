// BLS12-381 field prime
export const P =
    0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaabn
// BLS12-381 curve order
export const R = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001n

// n mod m
export function mod(n: bigint, m: bigint): bigint {
    return ((n % m) + m) % m
}

// Euclidean GCD extended binary algorithm
function gcd(u: bigint, v: bigint, x1: bigint, x2: bigint, p: bigint): bigint {
    if (!(u > 0n && v > 0n)) {
        throw new Error(`u,v=${u},${v} must be greater than 0`)
    }

    // Base cases
    if (u === 1n) return mod(x1, P)
    if (v === 1n) return mod(x2, P)
    if (u % 2n === 0n) {
        if (x1 % 2n === 0n) {
            return gcd(u >> 1n, v, x1 >> 1n, x2, p)
        } else {
            return gcd(u >> 1n, v, (x1 + p) >> 1n, x2, p)
        }
    }
    if (v % 2n === 0n) {
        if (x2 % 2n === 0n) {
            return gcd(u, v >> 1n, x1, x2 >> 1n, p)
        } else {
            return gcd(u, v >> 1n, x1, (x2 + P) >> 1n, p)
        }
    }
    if (u >= v) {
        return gcd(u - v, v, x1 - x2, x2, p)
    }
    if (u < v) {
        return gcd(u, v - u, x1, x2 - x1, p)
    }
    throw new Error(`GCD failed with u,v=${u},${v} x1,x2=${x1},${x2} p=${p}`)
}

export interface Field {
    equals(rhs: Field): boolean
    add(rhs: Field): Field
    sub(rhs: Field): Field
    mul(rhs: Field): Field
    toString(): string
    toJSON(): any
    neg(): Field
    inv(): Field
    // mulByScalar(c: bigint): Field
}

export class Fq implements Field {
    x: bigint

    constructor(x: bigint) {
        this.x = mod(x, P)
    }

    static fromTuple([x]: bigint[]): Fq {
        return new Fq(x)
    }

    static zero(): Fq {
        return new Fq(0n)
    }

    equals(rhs: Fq): boolean {
        return this.x === rhs.x
    }

    add(rhs: Fq): Fq {
        return new Fq(this.x + rhs.x)
    }

    sub(rhs: Fq): Fq {
        return new Fq(this.x - rhs.x)
    }

    neg(): Fq {
        return new Fq(-this.x)
    }

    mul(rhs: Fq): Fq {
        return new Fq(this.x * rhs.x)
    }

    mulByScalar(c: bigint): Fq {
        return new Fq(this.x * c)
    }

    mulByNonResidue(): Fq {
        return this
    }

    inv() {
        if (this.x === 0n) {
            throw new Error(`Inversion of zero`)
        }
        return new Fq(gcd(this.x, P, 1n, 0n, P))
    }

    toString() {
        return `(Fq ${this.x.toString()})`
    }

    toJSON() {
        return {
            x: this.x.toString(),
        }
    }
}

/// Quadratic extension to Fq
export class Fq2 implements Field {
    x: Fq
    y: Fq
    constructor(x: Fq, y: Fq) {
        this.x = x
        this.y = y
    }

    static fromTuple([u0, u1]: [bigint, bigint]): Fq2 {
        return new Fq2(new Fq(u0), new Fq(u1))
    }

    static fromNumber(n: bigint): Fq2 {
        return Fq2.fromTuple([n, 0n])
    }

    static zero(): Fq2 {
        return new Fq2(new Fq(0n), new Fq(0n))
    }

    static one(): Fq2 {
        return new Fq2(new Fq(1n), new Fq(0n))
    }

    equals(rhs: Fq2): boolean {
        return this.x.equals(rhs.x) && this.y.equals(rhs.y)
    }

    add(rhs: Fq2): Fq2 {
        return new Fq2(this.x.add(rhs.x), this.y.add(rhs.y))
    }

    sub(rhs: Fq2): Fq2 {
        return new Fq2(this.x.sub(rhs.x), this.y.sub(rhs.y))
    }

    neg(): Fq2 {
        return new Fq2(this.x.neg(), this.y.neg())
    }

    // TODO: Karatsuba algorithm
    mul(rhs: Fq2): Fq2 {
        const a0 = this.x
        const a1 = this.y
        const b0 = rhs.x
        const b1 = rhs.y
        const u = a0.mul(b0).sub(a1.mul(b1))
        const v = a1.mul(b0).add(a0.mul(b1))
        return new Fq2(u, v)
    }

    mulByScalar(c: bigint): Fq2 {
        const scalar = Fq2.fromNumber(c)
        return this.mul(scalar)
    }

    mulByNonResidue(): Fq2 {
        const x = this.x
        const y = this.y
        return new Fq2(x.sub(y), x.add(y))
    }

    inv(): Fq2 {
        const x0 = this.x
        const y0 = this.y
        const factor = x0.mul(x0).add(y0.mul(y0)).inv()
        const x1 = x0.mul(factor)
        const y1 = y0.neg().mul(factor)
        return new Fq2(x1, y1)
    }

    toString() {
        return `(Fq2 ${this.x.toString()}, ${this.y.toString()})`
    }

    toJSON() {
        return {
            x: this.x.toJSON(),
            y: this.y.toJSON(),
        }
    }

    toTuple(): [bigint, bigint] {
        return [this.x.x, this.y.x]
    }
}

/// Cubic extension to Fq2
export class Fq6 implements Field {
    x: Fq2
    y: Fq2
    z: Fq2
    constructor(x: Fq2, y: Fq2, z: Fq2) {
        this.x = x
        this.y = y
        this.z = z
    }

    static fromTuple([x, y, z]: [Fq2, Fq2, Fq2]): Fq6 {
        return new Fq6(x, y, z)
    }

    static fromNumber(n: bigint): Fq6 {
        return Fq6.fromTuple([Fq2.fromNumber(n), Fq2.zero(), Fq2.zero()])
    }

    static zero(): Fq6 {
        return Fq6.fromTuple([Fq2.zero(), Fq2.zero(), Fq2.zero()])
    }

    static one(): Fq6 {
        return Fq6.fromTuple([Fq2.one(), Fq2.zero(), Fq2.zero()])
    }

    equals(rhs: Fq6): boolean {
        return this.x.equals(rhs.x) && this.y.equals(rhs.y) && this.z.equals(rhs.z)
    }

    neg(): Fq6 {
        return Fq6.zero().sub(this)
    }

    add(rhs: Fq6): Fq6 {
        const a = this.x
        const b = this.y
        const c = this.z
        const x = rhs.x
        const y = rhs.y
        const z = rhs.z
        return new Fq6(a.add(x), b.add(y), c.add(z))
    }

    sub(rhs: Fq6): Fq6 {
        const a = this.x
        const b = this.y
        const c = this.z
        const x = rhs.x
        const y = rhs.y
        const z = rhs.z
        return new Fq6(a.sub(x), b.sub(y), c.sub(z))
    }

    // TODO: Optimise
    mul(rhs: Fq6): Fq6 {
        const a0 = this.x
        const a1 = this.y
        const a2 = this.z
        const b0 = rhs.x
        const b1 = rhs.y
        const b2 = rhs.z
        const t0 = a0.mul(b0)
        const t1 = a1.mul(b1)
        const t2 = a2.mul(b2)
        const z0 = t0.add(a1.add(a2).mul(b1.add(b2)).sub(t1.add(t2)).mulByNonResidue())
        const z1 = a0.add(a1).mul(b0.add(b1)).sub(t0.add(t1)).add(t2.mulByNonResidue())
        const z2 = t1.add(a0.add(a2).mul(b0.add(b2))).sub(t0.add(t2))
        return new Fq6(z0, z1, z2)
    }

    mulByNonResidue(): Fq6 {
        return Fq6.fromTuple([this.z.mulByNonResidue(), this.x, this.y])
    }

    inv(): Fq6 {
        const a0 = this.x
        const a1 = this.y
        const a2 = this.z
        const t0 = a0.mul(a0).sub(a1.mul(a2).mulByNonResidue())
        const t1 = a2.mul(a2).mulByNonResidue().sub(a0.mul(a1))
        const t2 = a1.mul(a1).sub(a0.mul(a2))
        const factor = a0
            .mul(t0)
            .add(a2.mul(t1).mulByNonResidue())
            .add(a1.mul(t2).mulByNonResidue())
            .inv()
        return Fq6.fromTuple([t0.mul(factor), t1.mul(factor), t2.mul(factor)])
    }

    toString() {
        return `(Fq6 ${this.x.toString()}, ${this.y.toString()}, ${this.z.toString()})`
    }

    toJSON() {
        return {
            x: this.x.toJSON(),
            y: this.y.toJSON(),
            z: this.z.toJSON(),
        }
    }

    toTuple(): [
        ReturnType<Fq2['toTuple']>,
        ReturnType<Fq2['toTuple']>,
        ReturnType<Fq2['toTuple']>,
    ] {
        return [this.x.toTuple(), this.y.toTuple(), this.z.toTuple()]
    }
}

/// 12th degree extension to Fq
export class Fq12 implements Field {
    x: Fq6
    y: Fq6
    constructor(x: Fq6, y: Fq6) {
        this.x = x
        this.y = y
    }

    static fromTuple([x, y]: [Fq6, Fq6]): Fq12 {
        return new Fq12(x, y)
    }

    static fromNumber(n: bigint): Fq12 {
        return Fq12.fromTuple([Fq6.fromNumber(n), Fq6.zero()])
    }

    static zero(): Fq12 {
        return Fq12.fromTuple([Fq6.zero(), Fq6.zero()])
    }

    static one(): Fq12 {
        return Fq12.fromTuple([Fq6.one(), Fq6.zero()])
    }

    equals(rhs: Fq12): boolean {
        return this.x.equals(rhs.x) && this.y.equals(rhs.y)
    }

    neg(): Fq12 {
        return Fq12.zero().sub(this)
    }

    add(rhs: Fq12): Fq12 {
        const a = this.x
        const b = this.y
        const x = rhs.x
        const y = rhs.y
        return new Fq12(a.add(x), b.add(y))
    }

    sub(rhs: Fq12): Fq12 {
        const a = this.x
        const b = this.y
        const x = rhs.x
        const y = rhs.y
        return new Fq12(a.sub(x), b.sub(y))
    }

    mul(rhs: Fq12): Fq12 {
        const a0 = this.x
        const a1 = this.y
        const b0 = rhs.x
        const b1 = rhs.y
        const t1 = a0.mul(b0)
        const t2 = a1.mul(b1)
        const x = t1.add(t2.mulByNonResidue())
        const y = a0.add(a1).mul(b0.add(b1)).sub(t1.add(t2))
        return new Fq12(x, y)
    }

    inv(): Fq12 {
        const a0 = this.x
        const a1 = this.y
        const factor = a0.mul(a0).sub(a1.mul(a1).mulByNonResidue()).inv()
        const x = a0.mul(factor)
        const y = a1.neg().mul(factor)
        return new Fq12(x, y)
    }

    toString() {
        return `(Fq12 ${this.x.toString()}, ${this.y.toString()})`
    }

    toJSON() {
        return {
            x: this.x.toJSON(),
            y: this.y.toJSON(),
        }
    }

    toTuple(): [ReturnType<Fq6['toTuple']>, ReturnType<Fq6['toTuple']>] {
        return [this.x.toTuple(), this.y.toTuple()]
    }
}
