- stress kernel calculation
- layered elastic models
   - Haskell matrix?

- point source model

- try `#pragma omp parallel for` (OpenMP) in okada_disloc3d.c and tde c files for the C functions

- use Sphinx to create the complete docs and upload to ReadTheDocs

- Github Actions to build and upload the packages, deployment failed 

- [2026-08-05 done] Linux 源码构建修复(readme.md/CMakelists.txt 大小写 + c0_ static),发布 0.1.3

- macOS 手动构建与发布步骤(下次发布新版本时,在 Mac 上操作):
  1. 升版本号:pyproject.toml 中 version = "x.y.z"
  2. 构建(注意:不要从仓库根目录运行,本地 build/ 目录会遮蔽 build 模块):
       cd /tmp && python -m build --wheel --sdist -o /tmp/dislocs_release <仓库路径>
  3. 预检:twine check /tmp/dislocs_release/*
  4. 上传:twine upload /tmp/dislocs_release/*
  5. 验证:新建 venv,pip install dislocs,跑冒烟测试(import dislocs + rde/tde)
  6. 提交:git commit -m "release: bump version to x.y.z"(不要带 Co-Authored-By)

- 可考虑:cibuildwheel 全平台 wheel(macOS arm64 + Linux x86_64,cp311/cp312/cp313),
  解决目前只有 macOS arm64/cp311 wheel、其他平台/版本需源码编译的问题

