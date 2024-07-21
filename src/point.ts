import { Field, Fq, Fq12, Fq2, Fq6, R, X } from './ff'
import { toBigInt } from './utils'

export interface Point<F extends Field> {
    x: F
    y: F
    equals(rhs: Point<F>): boolean
    isOnCurve(): boolean
    isInSubgroup(): boolean
    add(rhs: Point<F>): Point<F>
    double(): Point<F>
    neg(): Point<F>
    mul(c: bigint): Point<F>
    clone(): Point<F>
}

export interface PointClassType<F extends Field> {
    new (x: F, y: F): Point<F>
    generator(): Point<F>
    infinity(): Point<F>
}

export type PointInstanceType<F extends Field> = InstanceType<PointClassType<F>>

export class PointG1 implements PointInstanceType<Fq> {
    x: Fq
    y: Fq

    constructor(x: Fq, y: Fq) {
        this.x = x
        this.y = y
    }

    static fromBytes(bytes: Uint8Array): PointG1 {
        const firstByte = bytes[0]
        const isCompressed = ((firstByte >> 7) & 1) === 1
        const isInfinity = ((firstByte >> 6) & 1) === 1
        const yIsLexLargest = ((firstByte >> 5) & 1) === 1
        // Clear 3 MSB
        bytes[0] &= 0b00011111
        if (isCompressed && bytes.byteLength === 48) {
            // is compressed
            const x = new Fq(toBigInt(bytes))
            if (isInfinity) {
                if (!x.equals(Fq.zero())) {
                    throw new Error(`Invalid compressed point: expected zero, got ${x}`)
                }
                return PointG1.infinity()
            } else {
                // is not point at infinity; recover y from curve equation
                // y^2 = x^3 + 4 -> y = sqrt(x^3 + 4)
                let y = x.mul(x).mul(x).add(new Fq(4n)).sqrt()
                if (yIsLexLargest === y.signBigEndian()) {
                    y = y.neg()
                }
                return new PointG1(x, y)
            }
        } else if (!isCompressed && bytes.byteLength === 96) {
            // is uncompressed
            const x = new Fq(toBigInt(bytes.slice(0, 48)))
            const y = new Fq(toBigInt(bytes.slice(48, 96)))
            return new PointG1(x, y)
        } else {
            throw new Error(`Invalid point encoding for G1 (${bytes.byteLength} bytes)`)
        }
    }

    static generator(): PointG1 {
        // https://github.com/zcash/librustzcash/blob/6e0364cd42a2b3d2b958a54771ef51a8db79dd29/pairing/src/bls12_381/README.md#generators
        return new PointG1(
            new Fq(
                3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507n,
            ),
            new Fq(
                1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569n,
            ),
        )
    }

    static infinity(): PointG1 {
        return new PointG1(new Fq(0n), new Fq(0n))
    }

    equals(rhs: PointG1): boolean {
        return this.x.equals(rhs.x) && this.y.equals(rhs.y)
    }

    double(): PointG1 {
        // Base case: point at infinity
        if (this.equals(PointG1.infinity())) {
            return this
        }

        // NB: inv is expensive
        const x1 = this.x
        const y1 = this.y
        const k = x1.mul(x1).mulByScalar(3n).mul(y1.mulByScalar(2n).inv())
        const x = k.mul(k).sub(x1).sub(x1)
        const y = k.mul(x1.sub(x)).sub(y1)
        return new PointG1(x, y)
    }

    add(rhs: PointG1): PointG1 {
        // Base case: add identity
        if (this.equals(PointG1.infinity())) {
            return rhs
        }
        if (rhs.equals(PointG1.infinity())) {
            return this
        }

        const x1 = this.x
        const y1 = this.y
        const x2 = rhs.x
        const y2 = rhs.y
        if (x1.equals(x2) && y1.equals(y2)) {
            return this.double()
        } else if (x1.equals(x2) && !y1.equals(y2)) {
            return PointG1.infinity()
        } else {
            const k = y2.sub(y1).mul(x2.sub(x1).inv())
            const x = k.mul(k).sub(x1).sub(x2)
            const y = k.mul(x1.sub(x)).sub(y1)
            return new PointG1(x, y)
        }
    }

    neg(): PointG1 {
        return new PointG1(this.x, this.y.neg())
    }

    mul(c: bigint): PointG1 {
        // Base case: mul by zero
        if (c === 0n) return PointG1.infinity()
        // Case 1: invalid scalar
        // if (!(0n < c && c < R)) throw new Error(`Invalid scalar: ${c}`)
        // Case 2: mul by one
        if (c === 1n) return this.clone()

        // Double and add algorithm
        let base = this.clone()
        let acc = PointG1.infinity()
        while (c > 0n) {
            if (c % 2n === 1n) {
                // add remainder if not evenly divisible
                acc = acc.add(base)
            }
            base = base.double()
            c >>= 1n // Halve the scalar
        }
        return acc
    }

    clearCofactor(): PointG1 {
        // To map an element p in E(Fq) to G1, we can exponentiate p by 1-x
        // https://eprint.iacr.org/2019/403.pdf Section 5 "Clearing cofactors"
        return this.mul(X).add(this)
    }

    isOnCurve(): boolean {
        // Base case: point at infinity
        if (this.equals(PointG1.infinity())) {
            return false
        }

        // y^2 = x^3 + 4
        const lhs = this.y.mul(this.y)
        const rhs = this.x.mul(this.x).mul(this.x).add(new Fq(4n))
        return lhs.equals(rhs)
    }

    isInSubgroup(): boolean {
        return this.mul(R).equals(PointG1.infinity())
    }

    clone(): PointG1 {
        return new PointG1(this.x, this.y)
    }
}

export class PointG2 implements PointInstanceType<Fq2> {
    x: Fq2
    y: Fq2

    constructor(x: Fq2, y: Fq2) {
        this.x = x
        this.y = y
    }

    static fromBytes(bytes: Uint8Array): PointG2 {
        const firstByte = bytes[0]
        const isCompressed = ((firstByte >> 7) & 1) === 1
        const isInfinity = ((firstByte >> 6) & 1) === 1
        const yIsLexLargest = ((firstByte >> 5) & 1) === 1
        // Clear 3 MSB
        bytes[0] &= 0b00011111
        if (isCompressed && bytes.byteLength === 96) {
            // is compressed
            const x = new Fq(toBigInt(bytes))
            if (isInfinity) {
                if (!x.equals(Fq.zero())) {
                    throw new Error(`Invalid compressed point: expected zero, got ${x}`)
                }
                return PointG2.infinity()
            } else {
                const c1 = bytes.slice(0, 48)
                const c0 = bytes.slice(48, 96)
                const x = Fq2.fromTuple([toBigInt(c0), toBigInt(c1)])
                // is not point at infinity; recover y from curve equation
                // y^2 = x^3 + 4 -> y = sqrt(x^3 + 4)
                let y = x
                    .mul(x)
                    .mul(x)
                    .add(Fq2.fromTuple([4n, 4n]))
                    .sqrt()
                if (yIsLexLargest === y.signBigEndian()) {
                    y = y.neg()
                }
                return new PointG2(x, y)
            }
        } else if (!isCompressed && bytes.byteLength === 192) {
            // is uncompressed
            const x = Fq2.fromTuple([toBigInt(bytes.slice(0, 48)), toBigInt(bytes.slice(48, 96))])
            const y = Fq2.fromTuple([
                toBigInt(bytes.slice(96, 144)),
                toBigInt(bytes.slice(144, 192)),
            ])
            return new PointG2(x, y)
        } else {
            throw new Error(`Invalid point encoding for G1 (${bytes.byteLength} bytes)`)
        }
    }

    static generator(): PointG2 {
        // https://github.com/zcash/librustzcash/blob/6e0364cd42a2b3d2b958a54771ef51a8db79dd29/pairing/src/bls12_381/README.md#generators
        return new PointG2(
            Fq2.fromTuple([
                352701069587466618187139116011060144890029952792775240219908644239793785735715026873347600343865175952761926303160n,
                3059144344244213709971259814753781636986470325476647558659373206291635324768958432433509563104347017837885763365758n,
            ]),
            Fq2.fromTuple([
                1985150602287291935568054521177171638300868978215655730859378665066344726373823718423869104263333984641494340347905n,
                927553665492332455747201965776037880757740193453592970025027978793976877002675564980949289727957565575433344219582n,
            ]),
        )
    }

    static infinity(): PointG2 {
        return new PointG2(new Fq2(new Fq(0n), new Fq(0n)), new Fq2(new Fq(0n), new Fq(0n)))
    }

    equals(rhs: PointG2): boolean {
        return this.x.equals(rhs.x) && this.y.equals(rhs.y)
    }

    double(): PointG2 {
        // Base case: point at infinity
        if (this.equals(PointG2.infinity())) {
            return this
        }

        // NB: inv is expensive
        const x1 = this.x
        const y1 = this.y
        const k = x1.mul(x1).mulByScalar(3n).mul(y1.mulByScalar(2n).inv())
        const x = k.mul(k).sub(x1).sub(x1)
        const y = k.mul(x1.sub(x)).sub(y1)
        return new PointG2(x, y)
    }

    add(rhs: PointG2): PointG2 {
        // Base case: add identity
        if (this.equals(PointG2.infinity())) {
            return rhs
        }
        if (rhs.equals(PointG2.infinity())) {
            return this
        }

        const x1 = this.x
        const y1 = this.y
        const x2 = rhs.x
        const y2 = rhs.y
        if (x1.equals(x2) && y1.equals(y2)) {
            return this.double()
        } else if (x1.equals(x2) && !y1.equals(y2)) {
            return PointG2.infinity()
        } else {
            const k = y2.sub(y1).mul(x2.sub(x1).inv())
            const x = k.mul(k).sub(x1).sub(x2)
            const y = k.mul(x1.sub(x)).sub(y1)
            return new PointG2(x, y)
        }
    }

    neg(): PointG2 {
        return new PointG2(this.x, this.y.neg())
    }

    mul(c: bigint): PointG2 {
        // Base case: mul by zero
        if (c === 0n) return PointG2.infinity()
        // Case 1: invalid scalar
        // if (!(0n < c && c < R)) throw new Error(`Invalid scalar: ${c}`)
        // Case 2: mul by one
        if (c === 1n) return this.clone()

        // Double and add algorithm
        let base = this.clone()
        let acc = PointG2.infinity()
        while (c > 0n) {
            if (c % 2n === 1n) {
                // add remainder if not evenly divisible
                acc = acc.add(base)
            }
            base = base.double()
            c >>= 1n // Halve the scalar
        }
        return acc
    }

    /// Untwist an affine point in G2 defined over Fq2 to a point defined over Fq12
    untwist(): PointGT {
        const x1 = this.x
        const y1 = this.y
        const root = Fq6.fromTuple([Fq2.zero(), Fq2.one(), Fq2.zero()])
        // wideX = [ Fq12 (Fq6 x1 0 0) 0 ] * inv (Fq12 root 0)
        const wideX = Fq12.fromTuple([Fq6.fromTuple([x1, Fq2.zero(), Fq2.zero()]), Fq6.zero()]).mul(
            Fq12.fromTuple([root, Fq6.zero()]).inv(),
        )
        // wideY = [ Fq12 (Fq6 y1 0 0) 0 ] * inv (Fq12 0 root)
        const wideY = Fq12.fromTuple([Fq6.fromTuple([y1, Fq2.zero(), Fq2.zero()]), Fq6.zero()]).mul(
            Fq12.fromTuple([Fq6.zero(), root]).inv(),
        )
        return new PointGT(wideX, wideY)
    }

    clearCofactor(): PointG2 {
        // G2 cofactor
        const h2 =
            0x5d543a95414e7f1091d50792876a202cd91de4547085abaa68a205b2e5a7ddfa628f1cb4d9e82ef21537e293a6691ae1616ec6e786f0c70cf1c38e31c7238e5n
        return this.mul(h2)
    }

    isOnCurve(): boolean {
        // Base case: point at infinity
        if (this.equals(PointG2.infinity())) {
            return false
        }

        // y^2 = x^3 + 4
        const lhs = this.y.mul(this.y)
        const rhs = this.x
            .mul(this.x)
            .mul(this.x)
            .add(Fq2.fromTuple([4n, 4n])) // ??????
        return lhs.equals(rhs)
    }

    isInSubgroup(): boolean {
        return this.mul(R).equals(PointG2.infinity())
    }

    clone(): PointG2 {
        return new PointG2(this.x, this.y)
    }
}

/// Point in G12
export class PointGT implements PointInstanceType<Fq12> {
    x: Fq12
    y: Fq12

    constructor(x: Fq12, y: Fq12) {
        this.x = x
        this.y = y
    }

    static generator(): PointGT {
        throw new Error('Not implemented')
    }

    static infinity(): PointGT {
        return new PointGT(Fq12.zero(), Fq12.zero())
    }

    equals(rhs: PointGT): boolean {
        return this.x.equals(rhs.x) && this.y.equals(rhs.y)
    }

    isOnCurve(): boolean {
        // Base case: point at infinity
        if (this.equals(PointGT.infinity())) {
            return false
        }

        // y^2 = x^3 + 4
        const lhs = this.y.mul(this.y)
        const rhs = this.x.mul(this.x).mul(this.x).add(Fq12.fromNumber(4n))
        return lhs.equals(rhs)
    }

    isInSubgroup(): boolean {
        throw new Error('Not implemented')
    }

    add(rhs: PointGT): PointGT {
        throw new Error('Not implemented')
    }

    double(): PointGT {
        throw new Error('Not implemented')
    }

    neg(): PointGT {
        throw new Error('Not implemented')
    }

    mul(c: bigint): PointGT {
        throw new Error('Not implemented')
    }

    clone(): PointGT {
        return new PointGT(this.x, this.y)
    }
}
