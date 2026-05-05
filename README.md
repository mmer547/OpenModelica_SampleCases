# OpenModelica_SampleCases

OpenModelica で動作確認するためのサンプルケース集です。  
現時点では、冷媒圧縮機の簡易 1 室モデル `Compression_Chamber` を収録しています。

## 収録モデル

### `Compression_Chamber`

`Compression_Chamber/Compression_Chamber.mo` に、以下の構成が含まれます。

- `Components`
  - `CompressionChamber`: 圧縮室（質量・エネルギー収支、周期体積変化）
  - `SuctionChamber`: 吸入側境界、主吸入流量、漏れ流量
  - `DischargeValve`: 吐出弁（一方向流れの簡易モデル）
  - `DischargeBoundary`: 吐出圧力境界
  - `SystemParameters`: 共通パラメータ（圧力・温度・面積・流量係数など）
- `Media`
  - `CO2`, `R32`, `R410A`（`ExternalMedia.Media.CoolPropMedium` を利用）
- `Interfaces`
  - `port_in`, `port_out`, `port_leak` などのコネクタ定義
- `Examples`
  - `SimpleChamberTest`: 主要コンポーネントを接続したテストモデル

## 実行方法（OpenModelica）

1. OpenModelica（OMEdit）を起動
2. `Compression_Chamber/Compression_Chamber.mo` を読み込み
3. クラスブラウザから `Compression_Chamber.Examples.SimpleChamberTest` を開く
4. シミュレーション実行（必要に応じて stop time / solver を調整）
5. `compressionChamber.p`, `compressionChamber.T`, `dischargeValve.m_flow_out` などをプロット

## 依存関係

- `Modelica` 4.1.0
- `ExternalMedia`（CoolProp backend を利用）

`Media` パッケージは `ExternalMedia.Media.CoolPropMedium` を拡張しています。  
環境によっては、OpenModelica 側で ExternalMedia/CoolProp の導入と設定が必要です。

## 追加予定

このリポジトリはサンプルケースを追加していく前提です。  
現状は `Compression_Chamber` のみを収録しています。