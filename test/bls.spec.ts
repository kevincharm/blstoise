import { expect } from 'chai'
import { expandMsgXmd, hashToField } from '../src/bls'

describe('BLS', () => {
    const dst = new Uint8Array(Buffer.from('BLS_SIG_BLS12381G1_XMD:SHA-256_SSWU_RO_NUL_', 'ascii'))
    const messageBytes = Buffer.from(BigInt(9162609n).toString(16).padStart(16, '0'), 'hex')

    it('expandMsgXmd', () => {
        const expanded = expandMsgXmd(dst, messageBytes, 128)
        expect(expanded).to.deep.equal(
            Buffer.from(
                '64bd67fd58460b4bcf582bb7e4416c50047205f514683ffec11e8f0fbd2b87abeb087e50c21107dd2bdc72a3b99a83997dad227b7265addb8ea69391ce4464c262b141416688c0a348b37d4ef66c28485e11058517161603948f034d1201a57f95b612be6ba90c527fa545c9e0ab14cc5a48300bc2b18cde1712169053b70031',
                'hex',
            ),
        )
    })

    it('hashToField', () => {
        const els = hashToField(dst, messageBytes, 2)
        expect(els).to.deep.equal([
            3183467405108202357042404752704076250514358716694473968260566636730176549897642760246105978651082651830670222904862n,
            199803644312870517767392274353217167986279151073684446806214433239694194512637470382804262477277127014069748910675n,
        ])
    })

    // it('verify drand beacon', () => {
    //     const publicKeyBytes = Buffer.from(
    //         '83cf0f2896adee7eb8b5f01fcad3912212c437e0073e911fb90022d3e760183c8c4b450b6a0a6c3ac6a5776a2d1064510d1fec758c921cc22b0e17e63aaf4bcb5ed66304de9cf809bd274ca73bab4af5a6e9c76a4bc09e76eae8991ef5ece45a',
    //         'hex',
    //     )
    //     const signatureBytes = Buffer.from(
    //         '9819b31d0aedebe12a414ab4405dfce28d10fb33d74dfd8b921dd09a9915626faa7c35551584a7255dd1658062931104',
    //         'hex',
    //     )
    //     const messageBytes = BigInt(9162609n).toString(16).padStart(16, '0') // pad to 64 bits == 8 bytes
    // })
})
