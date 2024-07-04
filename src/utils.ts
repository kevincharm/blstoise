import { randomBytes } from 'node:crypto'
import { Fq, Fq12, Fq2, Fq6, P } from './ff'

export function randomFq12(): Fq12 {
    return Fq12.fromTuple([randomFq6(), randomFq6()])
}

export function randomFq6(): Fq6 {
    return Fq6.fromTuple([randomFq2(), randomFq2(), randomFq2()])
}

export function randomFq2(): Fq2 {
    return new Fq2(randomFq(), randomFq())
}

export function randomFq(): Fq {
    return new Fq(randomBigIntModP(48, P))
}

export function randomBigIntModP(bytes: number, p: bigint): bigint {
    const randUpperBound = 2n ** BigInt(bytes * 8)
    const upperBound = randUpperBound - (randUpperBound % p)
    let rand: bigint
    while ((rand = randomBigInt(bytes)) >= upperBound) {}
    return rand
}

export function randomBigInt(bytes: number): bigint {
    return BigInt(`0x${randomBytes(bytes).toString('hex')}`)
}
