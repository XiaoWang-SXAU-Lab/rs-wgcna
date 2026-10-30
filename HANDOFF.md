# rs-wgcna 引擎交接（2026-10-30）

本批修复必须与 `wgcna-dashboard` 仓库 `main` 的同批防御校验一起使用。
当前引擎分支是 **master**，Dashboard 分支是 **main**。

## 到家后必须重新编译

请先同步两个仓库，再在相邻的 `wgcna-dashboard` 中执行隔离构建：

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\build-windows.ps1
```

脚本会重新编译 Rust worker，再构建并内嵌前端与 worker 到 GUI。
**不能继续用旧 EXE/App，也不能只更新前端或只复制旧 worker。**
详细发布目录、R 环境保留及验收清单以 Dashboard 的 HANDOFF.md / BUILDING.md 为准。
不要在同步源码目录里直接产生 node_modules / target 等跨平台依赖。

## 本批引擎修改

- 新增 `--convert-workbook <input.xls> <output.tsv>`，供 GUI 预检旧版 Excel。
- 损坏/类型错误的 params.json 明确拒绝，不再静默切换为默认参数。
- 归一化产生 NaN/Inf 时，在写入分析文件前指出具体基因、样本。
- 过滤后不足 2 个可用基因时，明确提示筛选设置问题，避免进入聚类 panic。
- 样本信息 one-hot 拒绝空值、重复/缺失样本，以及保留样本中仅一个类别。
- 样本删除同步、原始输入备份、参数范围及依赖检查由新版 Dashboard 统一处理。

## 已验证

macOS 隔离构建通过；以下 Rust 回归测试各 3 项通过：

```bash
cargo test --locked --release --bin rs-wgcna-worker defensive_tests
cargo test --locked --release --test trait_alignment
```

保留原有 cutree 编译警告，本批未改变网络/模块检测数学算法。
Windows 的原生程序与文件对话框仍需在家中机器重新编译后实际验收。
