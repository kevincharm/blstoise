import { expect } from 'chai'
import fs from 'node:fs/promises'
import path from 'node:path'
import { Fq, Fq12, Fq2, Fq6, Fr } from '../src/ff'
import { randomFr } from '../src/utils'

async function readTestVectors<S, P>(
    name: string,
    mapFn: (value: S, index: number) => P,
): Promise<P[]> {
    const vectors = await fs.readFile(path.resolve(__dirname, 'vectors', `${name}.json`), {
        encoding: 'utf-8',
    })
    return JSON.parse(vectors).map(mapFn) as P[]
}

// Primitive types
type PFq = bigint
type PFq2 = [PFq, PFq]
type PFq6 = [PFq2, PFq2, PFq2]
type PFq12 = [PFq6, PFq6]
type PFqVector = [a: PFq, b: PFq, result: PFq]
type PFq2Vector = [a: PFq2, b: PFq2, result: PFq2]
type PFq6Vector = [a: PFq6, b: PFq6, result: PFq6]
type PFq12Vector = [a: PFq12, b: PFq12, result: PFq12]

// Serialised types
type SFq = string
type SFq2 = [SFq, SFq]
type SFq6 = [SFq2, SFq2, SFq2]
type SFq12 = [SFq6, SFq6]
type SFqVector = [a: SFq, b: SFq, result: SFq]
type SFq2Vector = [a: SFq2, b: SFq2, result: SFq2]
type SFq6Vector = [a: SFq6, b: SFq6, result: SFq6]
type SFq12Vector = [a: SFq12, b: SFq12, result: SFq12]

// Deserialisers
const reviveFq = ([x, y, z]: SFq[]): PFqVector => [BigInt(x), BigInt(y), BigInt(z)]
const reviveFq2 = ([x, y, z]: SFq2[]): PFq2Vector => [
    [BigInt(x[0]), BigInt(x[1])],
    [BigInt(y[0]), BigInt(y[1])],
    [BigInt(z[0]), BigInt(z[1])],
]
const reviveFq6 = ([x, y, z]: SFq6[]): PFq6Vector => [
    [
        [BigInt(x[0][0]), BigInt(x[0][1])],
        [BigInt(x[1][0]), BigInt(x[1][1])],
        [BigInt(x[2][0]), BigInt(x[2][1])],
    ],
    [
        [BigInt(y[0][0]), BigInt(y[0][1])],
        [BigInt(y[1][0]), BigInt(y[1][1])],
        [BigInt(y[2][0]), BigInt(y[2][1])],
    ],
    [
        [BigInt(z[0][0]), BigInt(z[0][1])],
        [BigInt(z[1][0]), BigInt(z[1][1])],
        [BigInt(z[2][0]), BigInt(z[2][1])],
    ],
]
const reviveFq12 = ([x, y, z]: SFq12[]): PFq12Vector => [
    [
        [
            [BigInt(x[0][0][0]), BigInt(x[0][0][1])],
            [BigInt(x[0][1][0]), BigInt(x[0][1][1])],
            [BigInt(x[0][2][0]), BigInt(x[0][2][1])],
        ],
        [
            [BigInt(x[1][0][0]), BigInt(x[1][0][1])],
            [BigInt(x[1][1][0]), BigInt(x[1][1][1])],
            [BigInt(x[1][2][0]), BigInt(x[1][2][1])],
        ],
    ],
    [
        [
            [BigInt(y[0][0][0]), BigInt(y[0][0][1])],
            [BigInt(y[0][1][0]), BigInt(y[0][1][1])],
            [BigInt(y[0][2][0]), BigInt(y[0][2][1])],
        ],
        [
            [BigInt(y[1][0][0]), BigInt(y[1][0][1])],
            [BigInt(y[1][1][0]), BigInt(y[1][1][1])],
            [BigInt(y[1][2][0]), BigInt(y[1][2][1])],
        ],
    ],
    [
        [
            [BigInt(z[0][0][0]), BigInt(z[0][0][1])],
            [BigInt(z[0][1][0]), BigInt(z[0][1][1])],
            [BigInt(z[0][2][0]), BigInt(z[0][2][1])],
        ],
        [
            [BigInt(z[1][0][0]), BigInt(z[1][0][1])],
            [BigInt(z[1][1][0]), BigInt(z[1][1][1])],
            [BigInt(z[1][2][0]), BigInt(z[1][2][1])],
        ],
    ],
]

describe('finite fields', () => {
    describe('Fr', () => {
        it('random', () => {
            const r = randomFr()
            expect(r.value).to.not.equal(0n)
        })

        it('inv', () => {
            const expected: [Fr, Fr][] = [
                [1n, 1n],
                [
                    2n,
                    26217937587563095239723870254092982918845276250263818911301829349969290592257n,
                ],
                [
                    3n,
                    34957250116750793652965160338790643891793701667018425215069105799959054123009n,
                ],
            ].map(([x, y]) => [new Fr(x), new Fr(y)])
            for (const [x, y] of expected) {
                expect(x.inv().equals(y)).to.equal(true)
            }
        })
    })

    describe('Fq', () => {
        it('add', async () => {
            const testVectors = await readTestVectors<SFqVector, PFqVector>('fp_add', reviveFq)
            for (const [x, y, z] of testVectors) {
                const a = new Fq(x)
                const b = new Fq(y)
                expect(a.add(b).value).to.equal(z)
            }
        })

        it('sub', async () => {
            const testVectors = await readTestVectors<SFqVector, PFqVector>('fp_sub', reviveFq)
            for (const [x, y, z] of testVectors) {
                const a = new Fq(x)
                const b = new Fq(y)
                expect(a.sub(b).value).to.equal(z)
            }
        })

        it('mul', async () => {
            const testVectors = await readTestVectors<SFqVector, PFqVector>('fp_mul', reviveFq)
            for (const [x, y, z] of testVectors) {
                const a = new Fq(x)
                const b = new Fq(y)
                expect(a.mul(b).value).to.equal(z)
            }
        })

        it('inv', async () => {
            const testVectors = await readTestVectors<[SFq, SFq], [PFq, PFq]>(
                'fp_inv',
                ([x, z]) => [BigInt(x), BigInt(z)],
            )
            for (const [x, z] of testVectors) {
                const a = new Fq(x)
                expect(a.inv().value).to.equal(z)
            }
        })
    })

    describe('Fq2', () => {
        it('add', async () => {
            const testVectors = await readTestVectors<SFq2Vector, PFq2Vector>('fp2_add', reviveFq2)
            for (const [x, y, z] of testVectors) {
                const a = Fq2.fromTuple(x)
                const b = Fq2.fromTuple(y)
                const result = a.add(b)
                expect(result.x.value).to.eq(z[0])
                expect(result.y.value).to.eq(z[1])
            }
        })

        it('sub', async () => {
            const testVectors = await readTestVectors<SFq2Vector, PFq2Vector>(
                'fp2_sub',
                ([x, y, z]) => [
                    [BigInt(x[0]), BigInt(x[1])],
                    [BigInt(y[0]), BigInt(y[1])],
                    [BigInt(z[0]), BigInt(z[1])],
                ],
            )
            for (const [x, y, z] of testVectors) {
                const a = Fq2.fromTuple(x)
                const b = Fq2.fromTuple(y)
                const result = a.sub(b)
                expect(result.x.value).to.eq(z[0])
                expect(result.y.value).to.eq(z[1])
            }
        })

        it('mul', async () => {
            const testVectors = await readTestVectors<SFq2Vector, PFq2Vector>('fp2_mul', reviveFq2)
            for (const [x, y, z] of testVectors) {
                const a = Fq2.fromTuple(x)
                const b = Fq2.fromTuple(y)
                const result = a.mul(b)
                expect(result.x.value).to.eq(z[0])
                expect(result.y.value).to.eq(z[1])
            }
        })
    })

    describe('Fq6', () => {
        it('add', async () => {
            const testVectors = await readTestVectors<SFq6Vector, PFq6Vector>('fp6_add', reviveFq6)
            for (const [x, y, z] of testVectors) {
                const a = Fq6.fromTuple([x[0], x[1], x[2]])
                const b = Fq6.fromTuple([y[0], y[1], y[2]])
                const result = a.add(b)
                expect(result.equals(Fq6.fromTuple([z[0], z[1], z[2]]))).to.eq(true)
            }
        })

        it('sub', async () => {
            const testVectors = await readTestVectors<SFq6Vector, PFq6Vector>('fp6_sub', reviveFq6)
            for (const [x, y, z] of testVectors) {
                const a = Fq6.fromTuple([x[0], x[1], x[2]])
                const b = Fq6.fromTuple([y[0], y[1], y[2]])
                const result = a.sub(b)
                expect(result.equals(Fq6.fromTuple([z[0], z[1], z[2]]))).to.eq(true)
            }
        })

        it('mul', async () => {
            const testVectors = await readTestVectors<SFq6Vector, PFq6Vector>('fp6_mul', reviveFq6)
            for (const [x, y, z] of testVectors) {
                const a = Fq6.fromTuple([x[0], x[1], x[2]])
                const b = Fq6.fromTuple([y[0], y[1], y[2]])
                const result = a.mul(b)
                expect(result.equals(Fq6.fromTuple([z[0], z[1], z[2]]))).to.eq(true)
            }
        })
    })

    describe('Fq12', () => {
        it('add', async () => {
            const testVectors = await readTestVectors<SFq12Vector, PFq12Vector>(
                'fp12_add',
                reviveFq12,
            )
            for (const [x, y, z] of testVectors) {
                const a = Fq12.fromTuple([
                    [x[0][0], x[0][1], x[0][2]],
                    [x[1][0], x[1][1], x[1][2]],
                ])
                const b = Fq12.fromTuple([
                    [y[0][0], y[0][1], y[0][2]],
                    [y[1][0], y[1][1], y[1][2]],
                ])
                const result = a.add(b)
                expect(
                    result.equals(
                        Fq12.fromTuple([
                            [z[0][0], z[0][1], z[0][2]],
                            [z[1][0], z[1][1], z[1][2]],
                        ]),
                    ),
                ).to.eq(true)
            }
        })

        it('sub', async () => {
            const testVectors = await readTestVectors<SFq12Vector, PFq12Vector>(
                'fp12_sub',
                reviveFq12,
            )
            for (const [x, y, z] of testVectors) {
                const a = Fq12.fromTuple([
                    [x[0][0], x[0][1], x[0][2]],
                    [x[1][0], x[1][1], x[1][2]],
                ])
                const b = Fq12.fromTuple([
                    [y[0][0], y[0][1], y[0][2]],
                    [y[1][0], y[1][1], y[1][2]],
                ])
                const result = a.sub(b)
                expect(
                    result.equals(
                        Fq12.fromTuple([
                            [z[0][0], z[0][1], z[0][2]],
                            [z[1][0], z[1][1], z[1][2]],
                        ]),
                    ),
                ).to.eq(true)
            }
        })

        it('mul', async () => {
            const testVectors = await readTestVectors<SFq12Vector, PFq12Vector>(
                'fp12_mul',
                reviveFq12,
            )
            for (const [x, y, z] of testVectors) {
                const a = Fq12.fromTuple([
                    [x[0][0], x[0][1], x[0][2]],
                    [x[1][0], x[1][1], x[1][2]],
                ])
                const b = Fq12.fromTuple([
                    [y[0][0], y[0][1], y[0][2]],
                    [y[1][0], y[1][1], y[1][2]],
                ])
                const result = a.mul(b)
                expect(
                    result.equals(
                        Fq12.fromTuple([
                            [z[0][0], z[0][1], z[0][2]],
                            [z[1][0], z[1][1], z[1][2]],
                        ]),
                    ),
                ).to.eq(true)
            }
        })
    })
})
