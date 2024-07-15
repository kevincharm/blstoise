/// BLS signature scheme
import { sha256 } from '@kevincharm/sharmander256'
import { validatePairing } from './pairing'
import { PointG1, PointG2 } from './point'
import { Fq, modp } from './ff'
import { toBigEndianBuffer, toBigInt } from './utils'
import { mapToPointSSWU } from './sswu'

/// Verify signatures on G1 from an already-hashed message
export function rawVerifyG1(
    publicKey: PointG2,
    signature: PointG1,
    hashedMessage: PointG1,
): boolean {
    const ps = [hashedMessage, signature]
    const qs = [publicKey.neg(), PointG2.generator()]
    return validatePairing(ps, qs)
}

export function hashToPoint(dst: Uint8Array, msg: Uint8Array): PointG1 {
    const [x, y] = hashToField(dst, msg, 2)
    const [p0, p1] = mapToPointSSWU(x)
    const p = new PointG1(new Fq(p0), new Fq(p1))
    const [q0, q1] = mapToPointSSWU(y)
    const q = new PointG1(new Fq(q0), new Fq(q1))
    const r = p.add(q).clearCofactor()
    return r
}

export function hashToField(dst: Uint8Array, msg: Uint8Array, count: number): bigint[] {
    const L = 64
    const _msg = expandMsgXmd(dst, msg, L * count)
    const els: bigint[] = []
    for (let i = 0; i < count; i++) {
        const el = modp(toBigInt(_msg.slice(i * L, i * L + L)))
        els.push(el)
    }
    return els
}

export function expandMsgXmd(dst: Uint8Array, message: Uint8Array, byteLength: number): Uint8Array {
    const ell = Math.ceil(byteLength / 32)
    const domainLen = dst.length
    if (byteLength > 65536 || ell > 255 || domainLen > 255) {
        throw new Error(`Invalid lengths: dst=${domainLen}, byteLength=${byteLength}, ell=${ell}`)
    }
    const zpad = new Uint8Array(64)
    const lenInBytes = new Uint8Array([(byteLength >> 8) & 0xff, byteLength & 0xff])
    const b_0 = new Uint8Array([...zpad, ...message, ...lenInBytes, 0, ...dst, domainLen])
    const b0 = sha256(b_0)
    let b_i = new Uint8Array([...b0, 1, ...dst, domainLen])
    let bi = sha256(b_i)
    const out = new Uint8Array(byteLength)
    for (let i = 1; i < ell; i++) {
        const xored = toBigEndianBuffer(toBigInt(b0) ^ toBigInt(bi), 32)
        b_i = new Uint8Array([...xored, 1 + i, ...dst, domainLen])
        out.set(bi, (i - 1) * 32)
        bi = sha256(b_i)
    }
    // final
    out.set(bi, (ell - 1) * 32)
    return out
}
