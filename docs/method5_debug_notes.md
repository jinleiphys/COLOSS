# Method 5 库仑散射调试笔记

## 1. 问题描述

Method 5 (Lagrange-Legendre mesh) 在**无库仑**情况下工作正常，但在**有库仑**情况下给出错误结果。

### 无库仑测试结果 (正确)

| 量 | Method 1 | Method 5 | 一致性 |
|---|---|---|---|
| f_born | (4.75, 0.63) | (4.75, 0.63) | ✓ |
| f_sc | (-4.63, -0.29) | (-4.63, -0.29) | ✓ |
| f_total | (0.117, 0.341) | (0.116, 0.339) | ✓ |
| S-matrix | (0.350, 0.223) | (0.354, 0.221) | ✓ |
| Direct vs Integral | - | 差异 ~7e-5 | ✓ |

### 有库仑测试结果 (修复前 - 错误)

| 量 | Method 1 | Method 5 |
|---|---|---|
| f_born | (6.12, 0.69) | (6.12, 0.69) ✓ |
| f_sc | (-6.05, 0.04) | (-2.88, -1.02) ✗ |
| f_total | (0.072, 0.726) | (3.24, -0.33) ✗ |

**关键观察**：
- f_born 完全正确，说明库仑波函数和短程势计算正确
- f_sc 错误，说明矩阵方程求解的 ψ_sc 有问题
- Direct match 和 Integral 方法差异巨大 (~3)，且随 Rmax 增大而增大

---

## 2. 理论公式 (来自论文)

### 波函数分解 (Eq. 11)
```
ψ_α(r) = e^{iσ_l} F_l(η, kr) + ψ_α^{sc}(r)
```

### 散射方程 (Eq. 13)
```
[E - H_α(r)] ψ_α^{sc}(r) = e^{iσ_l} Ṽ_N F_l(η, kr)
```

### 散射振幅 (Eq. 15)
```
f_α(k) = -2μ/(ℏ²k²) e^{-iσ_l} ∫ dr F_l(η, kr) Ṽ_N(r) ψ_α(r)
```

其中 Ṽ_N = V_N + V_C^S (短程势 = 核势 + 有限尺寸库仑 - 点库仑)

---

## 3. 尝试过的修改 (均失败)

### 尝试 1: 源项加 exp(iσ_l) 因子
- 根据 Eq. 13，源项应该有 e^{iσ_l} 因子
- **结果**: 无改善

### 尝试 2: f_sc 除以 exp(2iσ_l)
- 根据 Eq. 15 的相位因子
- **结果**: 无改善

### 尝试 3: Direct match 加相位因子
- 渐近形式 ψ_sc → k f e^{iσ_l} O_l^{(+)}
- **结果**: 无改善

### 尝试 4: 源项使用 V_short 并加相位
- **结果**: 仍然不对

---

## 4. 最终解决方案 (2025-12-01) ✓

### 关键发现

问题的根本原因是**相位因子的处理方式**。虽然论文公式 (Eq. 11-15) 明确写出了 e^{iσ_l} 因子，但在 Lagrange-Legendre 方法的具体实现中，这些因子会相互抵消。

### 正确的实现

#### 1. 矩阵 M 使用完整势能 V_full
```fortran
! V_full = V_nuc + V_coul_finite (包含有限尺寸库仑势)
vmod = vmod / coeff_kin  ! U = 2*mu*V/hbar^2
```

#### 2. 源项 b_vec 使用短程势 V_short，**不加** exp(iσ_l)
```fortran
! V_short = V_nuc + V_coul_finite - V_coul_point
! b = U_short * F_l * sqrt(lambda)
vmod_short = vmod_short / coeff_kin
b_vec(ir) = vmod_short * fc_loc(l) * sqrt(leg5_w(ir))
! 注意：不乘 exp(iu*cph(l))
```

#### 3. f_sc 计算**不除以** exp(2iσ_l)
```fortran
f_sc_int = -f_sc_int / ecm
! 注意：不除以 exp(2.d0*iu*cph(l))
```

#### 4. f_l_direct 计算**不除以** exp(2iσ_l)
```fortran
f_l_direct = phi_R / (k * hhat_R)
! 注意：不除以 exp(iu*cph(l))
```

#### 5. S 矩阵使用简单公式（与 Method 1 相同）
```fortran
! S = 1 + 2ik*f_l (不乘 exp(2iσ_l))
smat = 1.d0 + 2.d0 * iu * k * f_l
```

### 物理解释：为什么不需要 e^{-iσ_l} 因子

根据公式 (15)：
```
f_α(k) = -2μ/(ℏ²k²) e^{-iσ_l} ∫ dr F_l(η, kr) Ṽ_N(r) ψ_α(r)
```

积分外面应该有 `e^{-iσ_l}` 因子，但代码中没有加，结果却是对的。原因如下：

公式 (15) 中的 ψ_α 是**完整波函数**：
```
ψ_α = e^{iσ_l} F_l + ψ_sc
```

将其代入积分：
```
∫ F_l Ṽ_N ψ_α = ∫ F_l Ṽ_N (e^{iσ_l} F_l + ψ_sc)
              = e^{iσ_l} ∫ F_l Ṽ_N F_l + ∫ F_l Ṽ_N ψ_sc
```

乘以积分外的 `e^{-iσ_l}`：
- **Born 项**：`e^{-iσ_l} × e^{iσ_l} × ∫ F_l Ṽ_N F_l = ∫ F_l Ṽ_N F_l` （相位抵消！）
- **sc 项**：`e^{-iσ_l} × ∫ F_l Ṽ_N ψ_sc`

但在我们的实现中：
- 入射波选择为 `F_l`（不带 `e^{iσ_l}`）
- 因此散射波 `ψ_sc` 也不带相位因子
- 所以 `∫ F_l Ṽ_N ψ_sc` 本身就不需要 `e^{-iσ_l}` 修正

**结论**：我们计算的 `f_l` 实际上是 `e^{iσ_l} × f_α`，即**吸收了相位因子的散射振幅**。这与 S 矩阵公式 `S = 1 + 2ik×f_l`（不含 `e^{2iσ_l}`）完全一致。

这与 Method 1（Lagrange-Laguerre + 复数旋转）的处理方式一致。

---

## 5. 修复后的测试结果 ✓

### 有库仑测试 (L=0, p + Ca40, 20 MeV)

| 量 | Method 1 (ctheta=10) | Method 5 | 差异 |
|---|---|---|---|
| S-matrix (Re) | -0.386111 | -0.385529 | 0.0006 |
| S-matrix (Im) | 0.138203 | 0.134202 | 0.004 |
| f_l (Re) | 0.0724 | 0.0703 | 0.002 |
| f_l (Im) | 0.7262 | 0.7259 | 0.0003 |
| Direct vs Integral | - | ~0.0008 | ✓ |

### 无库仑测试

修复后无库仑情况仍然正确。

---

## 6. 代码修改位置

文件: `src/scatt_method5.f`

### 源项 (lines 359-364)
```fortran
! Source term: use SHORT-RANGE potential Ṽ_N = V_N + V_C^S
! b = U_short * F_l * sqrt(lambda)
vmod_short = vmod_short / coeff_kin  ! U_short = 2*mu*V_short/hbar^2
b_vec(ir) = vmod_short * fc_loc(l) * sqrt(leg5_w(ir))
```

### f_sc 计算 (line 530)
```fortran
f_sc_int = -f_sc_int / ecm
```

### f_l_direct (lines 475-480)
```fortran
f_l_direct = phi_R / (k * hhat_R)
```

### S 矩阵 (line 548)
```fortran
! S = 1 + 2ik*f_l (same formula as Method 1)
smat = 1.d0 + 2.d0 * iu * k * f_l
```

---

## 7. 总结

Method 5 库仑散射问题已解决。关键点：

1. **矩阵 M**: 使用完整势能 V_full = V_nuc + V_coul_finite
2. **源项 b**: 使用短程势 V_short = V_nuc + V_coul_finite - V_coul_point
3. **相位因子**: 不需要在源项、f_sc、f_l_direct 中添加 exp(iσ_l) 因子
4. **S 矩阵**: S = 1 + 2ik*f_l（与 Method 1 相同）

Method 5 现在可以正确处理有库仑和无库仑两种情况。
