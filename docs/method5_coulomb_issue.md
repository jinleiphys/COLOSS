# Method 5 库仑散射问题分析

## 1. 问题概述

Method 5 (Lagrange-Legendre mesh 方法) 在无库仑情况下工作正确，但在有库仑情况下给出错误结果。

## 2. 测试条件

- 入射粒子: 质子 (zp=1, massp=1)
- 靶核: Ca-40 (zt=20, masst=40)
- 入射能量: E_lab = 20 MeV
- 边界: R_max = 50 fm, nr = 150
- 势参数: 标准 Woods-Saxon 光学势

## 3. 散射振幅计算公式

总散射振幅分为两部分:
```
f_l = f_born + f_sc
```

其中:
- f_born = -1/E_cm × ∫ V_short × F_l² dr  (Born 近似项)
- f_sc = -1/E_cm × ∫ V_short × F_l × ψ_sc dr / exp(2iσ)  (散射修正项)

V_short = V_nuc + V_coul_finite - V_coul_point (短程势)

## 4. 数值结果对比

### 4.1 无库仑情况 (zp=0, zt=0)

| 量 | Method 1 | Method 5 | 一致性 |
|---|---|---|---|
| f_born | (4.747, 0.631) | (4.747, 0.631) | ✓ |
| f_sc | (-4.630, -0.290) | (-4.630, -0.290) | ✓ |
| f_total | (0.117, 0.341) | (0.117, 0.341) | ✓ |
| S-matrix | (0.350, 0.223) | (0.350, 0.223) | ✓ |

**结论: 无库仑时 Method 5 完全正确**

### 4.2 有库仑情况 (zp=1, zt=20)

| 量 | Method 1 | Method 5 | 一致性 |
|---|---|---|---|
| f_born | (6.121, 0.686) | (6.121, 0.686) | ✓ |
| f_sc | (-6.050, 0.037) | 待验证 | 待验证 |
| f_total | (0.072, 0.723) | 待验证 | 待验证 |
| S-matrix | (-0.381, 0.136) | 待验证 | 待验证 |

**结论: 有库仑时 f_born 正确，f_sc 需要验证**

## 5. 当前实现状态

### 5.1 矩阵中的势能

当前代码使用完整相互作用力 (V_full = V_nuc + V_coul_finite):
```fortran
! 在 solve_scatt_method5 中：
call compute_potential_at_r5(leg5_r(ir), para, ich, vmod)
! vmod 是 V_short = V_nuc + V_coul_finite - V_coul_point
! 加回 V_coul_point 得到 V_full
if (abs(z12) > 1.d-10) then
    vmod = vmod + e2 * z12 / leg5_r(ir)
endif
```

### 5.2 源项

源项使用 V_full × F_l × √λ:
```fortran
b_vec(ir) = vmod * fc_loc(l) * sqrt(leg5_w(ir))
```

### 5.3 散射振幅计算

f_born 和 f_sc 使用 V_short 计算 (与 Method 1 一致):
```fortran
! f_born: V_short * F_l^2 积分
call compute_potential_at_r5(leg5_r(ir), para, ich, vmod)
f_born = sum[V_short * F_l^2 * dr]

! f_sc: V_short * F_l * psi_sc 积分
f_sc = sum[V_short * F_l * psi_sc * dr]
```

## 6. 理论分析

### 6.1 Baye 方法的基本方程

对于有库仑的散射问题，Schrodinger 方程为:
```
[H_0 + V_full - E] ψ = 0
```

其中:
- H_0 = -ℏ²/(2μ) d²/dr² + ℏ²l(l+1)/(2μr²) (自由粒子动能)
- V_full = V_nuc + V_coul_finite (完整相互作用)

### 6.2 散射波方程

将波函数分解为 ψ = φ_inc + ψ_sc:
```
[H_0 + V_full - E] ψ_sc = -V_full × φ_inc
```

其中 φ_inc 是入射波 (库仑情况下是库仑波函数 F_l)。

### 6.3 边界条件

在 r = R 处，ψ_sc 满足出射波边界条件:
```
ψ'_sc(R) = γ_s × ψ_sc(R)
```

其中 γ_s = (G' + iF')/(G + iF) 是库仑出射波的对数导数。

## 7. Method 1 与 Method 5 的对比

### 7.1 Method 1 (Lagrange-Laguerre)

- 使用 Lagrange-Laguerre 基函数
- Hamiltonian 矩阵使用 V_full = V_nuc + V_coul_finite
- B_vec (源项) 使用 V_short × F_l × exp(iσ) × √λ
- f_sc 计算: -sum(X_vec × B_vec) / ecm / exp(2iσ)

### 7.2 Method 5 (Lagrange-Legendre)

- 使用 Lagrange-Legendre 基函数 (Baye 精确矩阵)
- 矩阵使用 V_full = V_nuc + V_coul_finite
- 源项使用 V_full × F_l × √λ
- f_sc 通过积分公式计算

## 8. 待解决问题

1. **相位因子**: Method 1 的 B_vec 包含 exp(iσ) 因子，Method 5 是否需要？
2. **f_sc 计算**: 积分公式是否需要除以 exp(2iσ)？
3. **势能选择**: 源项应该用 V_short 还是 V_full？

## 9. 相关代码文件

- `src/scatt_method5.f`: Method 5 主代码
- `src/scatt.f`: Method 1 代码 (solve_scatt 函数)
- `src/matrix_element.f`: Method 1 的 B_vec 计算 (cal_b 函数)
- `test/conv_m5_Coul_R50.in`: 有库仑测试输入
- `test/conv_m5_noCoul_R50.in`: 无库仑测试输入

## 10. 下一步

1. 运行测试验证当前实现
2. 检查 Method 1 的 exp(iσ) 因子是否需要在 Method 5 中添加
3. 确认 f_sc 计算公式中的相位因子处理
