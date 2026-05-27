# GrobCovBridge.jl

[English](README.md)

> **⚠ 実験的パッケージ — メンテナンス保証なし**
>
> このパッケージは一時的なブリッジ実装であり、Julia General レジストリには**登録されていません**。
> `grobcov.lib` の機能が OSCAR.jl に上流統合されるまでの暫定措置として存在します。
> 予告なく破壊的変更が行われる場合があります。プロダクションコードでの使用は推奨しません。

Singular の [`grobcov.lib`](https://www.singular.uni-kl.de/Manual/4-0-2/sing_956.htm) を
[OSCAR.jl](https://docs.oscar-system.org/) の多項式イデアルから呼び出すための Julia ブリッジパッケージです。

内部では `Singular.jl` の `low_level_caller_ring` を使って Singular インタプリタを呼び出し、
Oscar の `MPolyIdeal` と `grobcov.lib` が要求するパラメトリック環 `Q(params)[main_vars]`
の間の変換を自動的に処理します。

## インストール

```julia
using Pkg
Pkg.add(url="https://github.com/xiupos/GrobCovBridge.jl")
```

ローカルのクローンから開発モードでインストールする場合:

```julia
Pkg.develop(path="/path/to/GrobCovBridge.jl")
```

## クイックスタート

```julia
using Oscar, GrobCovBridge

R, (a, x, y) = polynomial_ring(QQ, ["a", "x", "y"])
I = ideal(R, [a*x - 1, a*y - 1])

# Groebner カバーを計算する。[a] がパラメータ変数、[x, y] が主変数。
result = grobcov(I, [a], [x, y])

# result は Vector{GrobCovSegment}
for seg in result
    println("lpp:     ", seg.lpp)
    println("basis:   ", gens(seg.basis))
    println("segment: ", seg.segment)
    println()
end
```

### 2 パラメータの例

```julia
R, (a, b, x, y) = polynomial_ring(QQ, ["a", "b", "x", "y"])
I = ideal(R, [a*x^2 + b*y, x - b*y])
result = grobcov(I, [a, b], [x, y])
```

### 包括的 Groebner システム（cgsdr）

```julia
R, (a, x) = polynomial_ring(QQ, ["a", "x"])
I = ideal(R, [a*x - 1])
result = cgsdr(I, [a], [x])

# result は Vector{CGSSegment}
for seg in result
    println("null:    ", gens(seg.null))
    println("nonnull: ", gens(seg.nonnull))
    println("basis:   ", gens(seg.basis))
    println()
end
```

## API

### `grobcov(I, params, main_vars; options=Int[])`

パラメトリックイデアル `I` の Groebner カバーを計算します。

- `I`: `QQ` 上の `Oscar.MPolyIdeal`
- `params`: `base_ring(I)` の生成元のうち、パラメータとして扱う変数のベクトル
- `main_vars`: 多項式変数として扱う生成元のベクトル
- `options`: Singular プロシージャに渡す整数フラグのリスト（省略可）

戻り値: `Vector{GrobCovSegment}`

#### `GrobCovSegment`

| フィールド | 型 | 説明 |
|-----------|-----|------|
| `lpp` | `Vector{MPolyRingElem}` | 基底の先頭冪積 |
| `basis` | `MPolyIdeal` | このセグメント上の被約 Groebner 基底 |
| `segment` | `Vector` | P 表現: `[[prime_ideal, [excl_ideal, ...]], ...]` |

`segment` はパラメータ空間の構成可能集合
`∪_i (V(prime_i) \ ∪_j V(excl_ij))` を表します。

### `cgsdr(I, params, main_vars; options=Int[])`

互いに素な被約表現による包括的 Groebner システムを計算します。

戻り値: `Vector{CGSSegment}`

#### `CGSSegment`

| フィールド | 型 | 説明 |
|-----------|-----|------|
| `null` | `MPolyIdeal` | セグメントのゼロ集合を定めるイデアル |
| `nonnull` | `MPolyIdeal` | 回避すべき集合を定めるイデアル |
| `basis` | `MPolyIdeal` | このセグメント上の Groebner 基底 |

## 動作要件

- Julia ≥ 1.10
- Oscar.jl ≥ 1.0（Singular.jl および `grobcov.lib` を含む）

## ライセンス

[GPL v3+](LICENSE)
