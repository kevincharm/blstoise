import { Fq, Fq12, Fq2, Fq6, Fr, Q, R } from './ff'

export function randomFq12(): Fq12 {
    return new Fq12(randomFq6(), randomFq6())
}

export function randomFq6(): Fq6 {
    return new Fq6(randomFq2(), randomFq2(), randomFq2())
}

export function randomFq2(): Fq2 {
    return new Fq2(randomFq(), randomFq())
}

export function randomFq(): Fq {
    return new Fq(randomBigIntModP(48, Q))
}

export function randomFr(): Fr {
    return new Fr(randomBigIntModP(31, R))
}

export function randomBigIntModP(bytes: number, p: bigint): bigint {
    const randUpperBound = 2n ** BigInt(bytes * 8)
    const upperBound = randUpperBound - (randUpperBound % p)
    let rand: bigint
    while ((rand = randomBigInt(bytes)) >= upperBound) {}
    return rand
}

export function randomBigInt(bytes: number): bigint {
    return toBigInt(crypto.getRandomValues(new Uint8Array(bytes)))
}

export function toBigInt(bytes: Uint8Array): bigint {
    return BigInt(`0x${toHex(bytes)}`)
}

export function toBigEndianBuffer(big: bigint, byteLength: number): Uint8Array {
    const buf = new Uint8Array(byteLength)
    let i = byteLength
    while (big > 0n) {
        buf[--i] = Number(big & 0xffn)
        big >>= 8n
    }
    return buf
}

export function toHex(bytes: Uint8Array) {
    return Array.from(bytes)
        .map((x) => x.toString(16).padStart(2, '0'))
        .join('')
}

export function assert(cond: boolean, msg?: string) {
    if (!cond) {
        if (msg) {
            throw new Error(msg)
        } else {
            throw new Error('Assertion failed')
        }
    }
}
