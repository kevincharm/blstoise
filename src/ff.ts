// BLS12-381 field prime
export const P =
    0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaabn
// BLS12-381 curve order
export const R = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001n
// BLS12-381 x parameter used for the construction of the curve
export const X = 0xd201000000010000n

// n mod m
export function mod(n: bigint, m: bigint): bigint {
    return ((n % m) + m) % m
}

export function modp(n: bigint): bigint {
    return mod(n, P)
}

export function egcd(a: bigint, b: bigint): [bigint, bigint, bigint] {
    if (a === 0n) {
        return [b, 0n, 1n]
    } else {
        const [g, y, x] = egcd(mod(b, a), a)
        return [g, x - (b / a) * y, y]
    }
}

// x * y (mod p)
export function mulmodp(x: bigint, y: bigint): bigint {
    return mod(x * y, P)
}

export function modexp(x: bigint, y: bigint, p: bigint): bigint {
    let result = 1n
    while (y > 0) {
        const lsb = y & 1n
        y = y / 2n
        if (lsb) {
            result = mod(result * x, p)
        }
        x = mod(x * x, p)
    }
    return result
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

function frobeniusCoeffsPowers(modulus: bigint, degree: bigint, num?: bigint, divisor?: bigint) {
    divisor = divisor || degree
    num = num || 1n
    const tower_modulus = modulus ** degree
    const powers = []
    for (let i = 0n; i < num; i++) {
        const a = i + 1n
        let q_power = 1n
        const ps = []
        for (let j = 0; j < degree; j++) {
            ps.push(mod((a * q_power - a) / divisor, tower_modulus))
            q_power *= modulus
        }
        powers.push(ps)
    }
    return powers
}

// https://gist.github.com/HarryR/eb5ad0e5de51633678e015a6b06969a1
function frobeniusCoeffs<F extends Field>(
    nonResidue: F,
    modulus: bigint,
    degree: bigint,
    num?: bigint,
    divisor?: bigint,
) {
    const coeffs: F[][] = []
    const allPowers = frobeniusCoeffsPowers(modulus, degree, num, divisor)
    for (let i = 0; i < allPowers.length; i++) {
        const powers = allPowers[i]
        coeffs.push([])
        for (const p_i of powers) {
            coeffs[i].push(nonResidue.exp(p_i) as F)
        }
    }
    return coeffs
}

export interface Field {
    equals(rhs: ThisType<this>): boolean
    add(rhs: ThisType<this>): ThisType<this>
    sub(rhs: ThisType<this>): ThisType<this>
    mul(rhs: ThisType<this>): ThisType<this>
    exp(rhs: bigint): ThisType<this>
    toString(): string
    toJSON(): any
    neg(): ThisType<this>
    inv(): ThisType<this>
    conjugate?(): ThisType<this>
    // mulByScalar(c: bigint):ThisType<this>
}

export class Fq implements Field {
    value: bigint

    constructor(x: bigint) {
        this.value = mod(x, P)
    }

    static zero(): Fq {
        return new Fq(0n)
    }

    equals(rhs: Fq): boolean {
        return this.value === rhs.value
    }

    lt(rhs: Fq): boolean {
        return this.value < rhs.value
    }

    add(rhs: Fq): Fq {
        return new Fq(this.value + rhs.value)
    }

    sub(rhs: Fq): Fq {
        return new Fq(this.value - rhs.value)
    }

    neg(): Fq {
        return new Fq(-this.value)
    }

    mul(rhs: Fq): Fq {
        return new Fq(this.value * rhs.value)
    }

    mulByScalar(c: bigint): Fq {
        return new Fq(this.value * c)
    }

    mulByNonResidue(): Fq {
        return this
    }

    inv() {
        if (this.value === 0n) {
            throw new Error(`Inversion of zero`)
        }
        return new Fq(gcd(this.value, P, 1n, 0n, P))
    }

    exp(y: bigint): Fq {
        let x: Fq = this
        let result = new Fq(1n)
        while (y > 0) {
            const lsb = y & 1n
            y = y / 2n
            if (lsb) {
                result = result.mul(x)
            }
            x = x.mul(x)
        }
        return result
    }

    legendre(): number {
        const x = this.exp((P - 1n) / 2n)
        if (x.equals(new Fq(P - 1n))) {
            return -1
        }
        if (!x.equals(Fq.zero()) && !x.equals(new Fq(1n))) {
            throw new Error(`Legendre failed: ${this}^{(${P}-1)//2} = ${x}`)
        }
        return Number(x)
    }

    sqrt(): Fq {
        const root = this.exp((P + 1n) / 4n)
        if (!root.mul(root).equals(this)) {
            throw new Error(`No square root exists for ${this}`)
        }
        return root
    }

    conjugate(): Fq {
        return new Fq(this.value)
    }

    signBigEndian(): boolean {
        const negV = this.neg()
        return this.lt(negV)
    }

    toString() {
        return `(Fq ${this.value.toString()})`
    }

    toJSON() {
        return {
            x: this.value.toString(),
        }
    }
}

/// Quadratic extension to Fq
export class Fq2 implements Field {
    // This is taken from consensus specs
    static readonly ORDER = P ** 2n
    static readonly EIGHTH_ROOTS_OF_UNITY = Array.from({ length: 8 }, (_, k) =>
        Fq2.fromTuple([1n, 1n]).exp((Fq2.ORDER * BigInt(k)) / 8n),
    )
    static readonly EVEN_EIGHT_ROOTS_OF_UNITY = Fq2.EIGHTH_ROOTS_OF_UNITY.filter(
        (_, i) => i % 2 === 0,
    )
    static readonly NON_RESIDUE = Fq2.fromTuple([1n, 1n])
    static readonly FROBENIUS_COEFFICIENTS = frobeniusCoeffs(new Fq(-1n), P, 2n)

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

    gt(rhs: Fq2): boolean {
        return this.x.value > rhs.x.value && this.y.value > rhs.y.value
    }

    lt(rhs: Fq2): boolean {
        return this.x.value < rhs.x.value && this.y.value < rhs.y.value
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

    exp(y: bigint): Fq2 {
        let x: Fq2 = this
        let result = Fq2.one()
        while (y > 0) {
            const lsb = y & 1n
            y = y / 2n
            if (lsb) {
                result = result.mul(x)
            }
            x = x.mul(x)
        }
        return result
    }

    // https://github.com/ethereum/trinity/blob/a1b0f058e7bc8e385c8dac3164c49098967fd5bb/eth2/_utils/bls.py#L63-L77
    sqrt(): Fq2 {
        const value = this
        const candidateSquareroot = value.exp((Fq2.ORDER + 8n) / 16n)
        const invValue = value.inv()
        const check = candidateSquareroot.mul(candidateSquareroot).mul(invValue)
        const rootIdx = Fq2.EVEN_EIGHT_ROOTS_OF_UNITY.findIndex((root) => root.equals(check))
        if (rootIdx === -1) throw new Error(`No square root exists for ${value}`)
        const root = Fq2.EIGHTH_ROOTS_OF_UNITY[rootIdx]
        const x1 = candidateSquareroot.mul(root.inv())
        const x2 = x1.neg()
        // (x1.coeffs[1], x1.coeffs[0]) > (x2.coeffs[1], x2.coeffs[0])
        if (x1.gt(x2)) {
            return x1
        } else {
            return x2
        }
    }

    conjugate(): Fq2 {
        return new Fq2(this.x, this.y.neg())
    }

    frobeniusMap(power: bigint): Fq2 {
        if (power & 1n) {
            return this.conjugate()
        }
        return new Fq2(this.x, this.y)
    }

    signBigEndian(): boolean {
        const negV = this.neg()
        return this.lt(negV)
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
        return [this.x.value, this.y.value]
    }
}

type Fq2Tuple = ReturnType<Fq2['toTuple']>

/// Cubic extension to Fq2
export class Fq6 implements Field {
    static readonly FROBENIUS_COEFFICIENTS = frobeniusCoeffs(Fq2.NON_RESIDUE, P, 6n, 2n, 3n)

    x: Fq2
    y: Fq2
    z: Fq2
    constructor(x: Fq2, y: Fq2, z: Fq2) {
        this.x = x
        this.y = y
        this.z = z
    }

    static fromTuple([x, y, z]: [Fq2Tuple, Fq2Tuple, Fq2Tuple]): Fq6 {
        return new Fq6(Fq2.fromTuple(x), Fq2.fromTuple(y), Fq2.fromTuple(z))
    }

    static fromNumber(n: bigint): Fq6 {
        return new Fq6(Fq2.fromNumber(n), Fq2.zero(), Fq2.zero())
    }

    static zero(): Fq6 {
        return new Fq6(Fq2.zero(), Fq2.zero(), Fq2.zero())
    }

    static one(): Fq6 {
        return new Fq6(Fq2.one(), Fq2.zero(), Fq2.zero())
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

    exp(y: bigint): Fq6 {
        let x: Fq6 = this
        let result = Fq6.one()
        while (y > 0) {
            const lsb = y & 1n
            y = y / 2n
            if (lsb) {
                result = result.mul(x)
            }
            x = x.mul(x)
        }
        return result
    }

    mulByNonResidue(): Fq6 {
        return new Fq6(this.z.mulByNonResidue(), this.x, this.y)
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
        return new Fq6(t0.mul(factor), t1.mul(factor), t2.mul(factor))
    }

    conjugate(): Fq6 {
        return new Fq6(this.x, this.y.neg(), this.z) // ???
    }

    frobeniusMap(power: bigint): Fq6 {
        const x = this.x.frobeniusMap(power)
        const c = Number(power % 6n)
        const y = this.y.frobeniusMap(power).mul(Fq6.FROBENIUS_COEFFICIENTS[0][c])
        const z = this.z.frobeniusMap(power).mul(Fq6.FROBENIUS_COEFFICIENTS[1][c])
        return new Fq6(x, y, z)
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

type Fq6Tuple = ReturnType<Fq6['toTuple']>

/// 12th degree extension to Fq
export class Fq12 implements Field {
    static readonly FROBENIUS_COEFFICIENTS = frobeniusCoeffs(Fq2.NON_RESIDUE, P, 12n, 1n, 6n)

    x: Fq6
    y: Fq6
    constructor(x: Fq6, y: Fq6) {
        this.x = x
        this.y = y
    }

    static fromTuple([x, y]: [Fq6Tuple, Fq6Tuple]): Fq12 {
        return new Fq12(Fq6.fromTuple(x), Fq6.fromTuple(y))
    }

    static fromNumber(n: bigint): Fq12 {
        return new Fq12(Fq6.fromNumber(n), Fq6.zero())
    }

    static zero(): Fq12 {
        return new Fq12(Fq6.zero(), Fq6.zero())
    }

    static one(): Fq12 {
        return new Fq12(Fq6.one(), Fq6.zero())
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

    exp(y: bigint): Fq12 {
        let x: Fq12 = this
        let result = Fq12.one()
        while (y > 0) {
            const lsb = y & 1n
            y = y / 2n
            if (lsb) {
                result = result.mul(x)
            }
            x = x.mul(x)
        }
        return result
    }

    inv(): Fq12 {
        const a0 = this.x
        const a1 = this.y
        const factor = a0.mul(a0).sub(a1.mul(a1).mulByNonResidue()).inv()
        const x = a0.mul(factor)
        const y = a1.neg().mul(factor)
        return new Fq12(x, y)
    }

    conjugate(): Fq12 {
        return new Fq12(this.x, this.y.neg())
    }

    frobeniusMap(power: bigint): Fq12 {
        const x = this.x.frobeniusMap(power)
        const { x: y0, y: y1, z: y2 } = this.y.frobeniusMap(power)
        const coeff = Fq12.FROBENIUS_COEFFICIENTS[0][Number(power % 12n)]
        return new Fq12(x, new Fq6(y0.mul(coeff), y1.mul(coeff), y2.mul(coeff)))
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

    toTuple(): [Fq6Tuple, Fq6Tuple] {
        return [this.x.toTuple(), this.y.toTuple()]
    }
}
