import os,sys
import numpy as np
from typing import Tuple, List, Optional
from scipy.spatial.distance import pdist, squareform
from itertools import combinations

class MolecularSurfaceCalculator:
    """分子体积和表面积计算器"""
    
    def __init__(self, molecule,method="grid",grid_spacing: float = 0.1):
        """
        初始化计算器
        
        参数:
            molecule: 分子对象
            method: grid or approximate
        """
        self.molecule = molecule
        self.coordinates = np.asarray(np.array(molecule.coordinates))
        self.elements = np.asarray(np.array(molecule.elements))
        self.n_atoms = len(self.elements)
        
        # 原子质量 (用于质心计算)
        #self.atomic_masses = self._get_atomic_masses()
        self.atomic_masses = np.array([atom.mass for atom in molecule.Atoms])
        
        # 原子范德华半径 (Å)
        #self.vdw_radii = self._get_vdw_radii()
        self.vdw_radii = np.array([atom.vdw_radius for atom in molecule.Atoms])
        
        # 计算质心
        self.center_of_mass = self._calculate_center_of_mass()
        
        # 将坐标平移到质心
        self.coords_centered = self.coordinates - self.center_of_mass
        self.method = method
        self.grid_spacing = grid_spacing
    
    def _calculate_center_of_mass(self) -> np.ndarray:
        """计算质心"""
        total_mass = np.sum(self.atomic_masses)
        com = np.sum(self.atomic_masses[:, np.newaxis] * self.coordinates, axis=0) / total_mass
        return com
    
    def _sphere_overlap_volume(self, r1: float, r2: float, d: float) -> float:
        """
        计算两个重叠球体的重叠体积
        
        参数:
            r1, r2: 两个球的半径
            d: 球心距离
        
        返回:
            重叠体积
        """
        if d >= r1 + r2:
            return 0.0
        if d <= abs(r1 - r2):
            return (4/3) * np.pi * min(r1, r2)**3
        
        R = r1
        r = r2
        
        part1 = (R**2 - r**2 + d**2) / (2 * d)
        part2 = (r**2 - R**2 + d**2) / (2 * d)
        
        volume = (np.pi * (R - part1)**2 * (d + part1 + r) / 3 +
                 np.pi * (r - part2)**2 * (d + part2 + R) / 3)
        
        return volume
    
    def _sphere_overlap_surface(self, r1: float, r2: float, d: float) -> float:
        """
        计算两个重叠球体的重叠表面积
        
        参数:
            r1, r2: 两个球的半径
            d: 球心距离
        
        返回:
            重叠表面积
        """
        if d >= r1 + r2:
            return 0.0
        if d <= abs(r1 - r2):
            return 4 * np.pi * min(r1, r2)**2
        
        h1 = r1 - (r1**2 - r2**2 + d**2) / (2 * d)
        h2 = r2 - (r2**2 - r1**2 + d**2) / (2 * d)
        
        overlap_area = 2 * np.pi * r1 * h1 + 2 * np.pi * r2 * h2
        return overlap_area
    
    def calculate_volume_approximate(self) -> float:
        """
        计算分子体积（近似方法，考虑原子间重叠）
        
        使用容斥原理近似计算
        """
        # 计算所有原子对的距离
        distances = squareform(pdist(self.coords_centered))
        
        # 计算单个球体体积之和
        single_sphere_volumes = (4/3) * np.pi * self.vdw_radii**3
        total_volume = np.sum(single_sphere_volumes)
        
        # 减去两两重叠的体积
        for i, j in combinations(range(self.n_atoms), 2):
            d = distances[i, j]
            overlap = self._sphere_overlap_volume(self.vdw_radii[i], self.vdw_radii[j], d)
            total_volume -= overlap
        
        return total_volume
    
    def calculate_surface_area_approximate(self) -> float:
        """
        计算分子表面积（近似方法，考虑原子间重叠）
        
        计算可访问表面积（Connolly surface的近似）
        """
        distances = squareform(pdist(self.coords_centered))
        
        total_area = 0.0
        
        for i in range(self.n_atoms):
            r_i = self.vdw_radii[i]
            accessible_fraction = 1.0
            
            # 计算每个原子被其他原子遮挡的比例
            for j in range(self.n_atoms):
                if i == j:
                    continue
                
                r_j = self.vdw_radii[j]
                d = distances[i, j]
                
                if d < r_i + r_j:
                    # 计算遮挡角度
                    if d <= abs(r_i - r_j):
                        # 完全被遮挡
                        if r_i > r_j:
                            accessible_fraction = 0.0
                            break
                    else:
                        # 部分遮挡
                        h_i = r_i - (r_i**2 - r_j**2 + d**2) / (2 * d)
                        cap_area = 2 * np.pi * r_i * h_i
                        accessible_fraction -= cap_area / (4 * np.pi * r_i**2)
            
            if accessible_fraction > 0:
                total_area += accessible_fraction * 4 * np.pi * r_i**2
        
        return max(total_area, 0.0)
    
    def calculate_volume_grid(self, grid_spacing: float = 0.1) -> Tuple[float, np.ndarray]:
        """
        使用网格方法计算分子体积
        
        参数:
            grid_spacing: 网格间距 (Å)
        
        返回:
            体积 (Å³) 和网格点数量
        """
        # 计算边界框
        coords = self.coords_centered
        radii = self.vdw_radii
        
        min_coord = np.min(coords - radii[:, np.newaxis], axis=0) - 0.5
        max_coord = np.max(coords + radii[:, np.newaxis], axis=0) + 0.5
        
        # 创建网格
        nx = int(np.ceil((max_coord[0] - min_coord[0]) / grid_spacing))
        ny = int(np.ceil((max_coord[1] - min_coord[1]) / grid_spacing))
        nz = int(np.ceil((max_coord[2] - min_coord[2]) / grid_spacing))
        
        x = np.linspace(min_coord[0], max_coord[0], nx)
        y = np.linspace(min_coord[1], max_coord[1], ny)
        z = np.linspace(min_coord[2], max_coord[2], nz)
        
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
        grid_points = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
        
        # 检查每个网格点是否在分子内
        inside = np.zeros(len(grid_points), dtype=bool)
        
        for i in range(self.n_atoms):
            dist = np.linalg.norm(grid_points - coords[i], axis=1)
            inside |= (dist <= radii[i])
        
        # 计算体积
        volume = np.sum(inside) * grid_spacing**3
        
        return volume, inside.reshape((nx, ny, nz))
    
    def calculate_surface_area_grid(self, grid_spacing: float = 0.1) -> float:
        """
        使用网格方法计算分子表面积
        
        参数:
            grid_spacing: 网格间距 (Å)
        
        返回:
            表面积 (Å²)
        """
        # 计算边界框
        coords = self.coords_centered
        radii = self.vdw_radii
        
        min_coord = np.min(coords - radii[:, np.newaxis], axis=0) - 0.5
        max_coord = np.max(coords + radii[:, np.newaxis], axis=0) + 0.5
        
        # 创建网格
        nx = int(np.ceil((max_coord[0] - min_coord[0]) / grid_spacing))
        ny = int(np.ceil((max_coord[1] - min_coord[1]) / grid_spacing))
        nz = int(np.ceil((max_coord[2] - min_coord[2]) / grid_spacing))
        
        x = np.linspace(min_coord[0], max_coord[0], nx)
        y = np.linspace(min_coord[1], max_coord[1], ny)
        z = np.linspace(min_coord[2], max_coord[2], nz)
        
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
        grid_points = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
        
        # 检查每个网格点是否在分子内
        inside = np.zeros(len(grid_points), dtype=bool)
        
        for i in range(self.n_atoms):
            dist = np.linalg.norm(grid_points - coords[i], axis=1)
            inside |= (dist <= radii[i])
        
        inside_grid = inside.reshape((nx, ny, nz))
        
        # 计算表面：找到至少有一个邻居在外部的内部点
        surface_area = 0.0
        surface_area_per_face = grid_spacing**2
        
        for i, j, k in zip(*np.where(inside_grid)):
            # 检查6个方向的邻居
            for di, dj, dk in [(-1,0,0), (1,0,0), (0,-1,0), (0,1,0), (0,0,-1), (0,0,1)]:
                ni, nj, nk = i + di, j + dj, k + dk
                if (ni < 0 or ni >= nx or nj < 0 or nj >= ny or 
                    nk < 0 or nk >= nz or not inside_grid[ni, nj, nk]):
                    surface_area += surface_area_per_face
        
        return surface_area

    def run(self,) -> Tuple[float, float]:
        """
        计算分子的体积和表面积
    
        参数:
            coordinates: 原子坐标数组, shape (n_atoms, 3), 单位: Angstrom
            elements: 原子元素符号列表
            method: 计算方法 ('grid' 或 'approximate')
            grid_spacing: 网格间距 (Å), 仅用于 'grid' 方法
    
        返回:
            volume: 分子体积 (Å³)
            surface_area: 分子表面积 (Å²)
        """
        #calculator = MolecularSurfaceCalculator(coordinates, elements)
    
        if self.method == 'grid':
            volume, _ = self.calculate_volume_grid(self.grid_spacing)
            surface_area = self.calculate_surface_area_grid(self.grid_spacing)
        elif self.method == 'approximate':
            volume = self.calculate_volume_approximate()
            surface_area = self.calculate_surface_area_approximate()
        else:
            raise ValueError(f"未知方法: {self.method}")
        
        self.molecule.volume = volume
        self.molecule.surface = surface_area
    
        return self.molecule


def calculate_moment_of_inertia(molecule,) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    计算分子的转动惯量张量和主转动惯量
    
    参数:
        coordinates: 原子坐标数组, shape (n_atoms, 3), 单位: Angstrom
        masses: 原子质量数组, shape (n_atoms,), 单位: atomic mass units (amu)
    
    返回:
        inertia_tensor: 转动惯量张量, shape (3, 3), 单位: amu * Å²
        principal_moments: 主转动惯量, shape (3,), 单位: amu * Å²
        eigenvectors: 主轴方向向量, shape (3, 3)
    """
    coordinates = np.asarray(np.array(molecule.coordinates))
    masses = np.asarray(np.array([atom.mass for atom in molecule.Atoms]))
    
    if coordinates.shape[0] != masses.shape[0]:
        raise ValueError("坐标和质量数组的长度不匹配")
    
    n_atoms = len(masses)
    
    # 计算质心
    total_mass = np.sum(masses)
    center_of_mass = np.sum(masses[:, np.newaxis] * coordinates, axis=0) / total_mass
    
    # 将坐标平移到质心
    coords_centered = coordinates - center_of_mass
    
    # 计算转动惯量张量
    I = np.zeros((3, 3))
    
    for i in range(n_atoms):
        m = masses[i]
        x, y, z = coords_centered[i]
        
        r_sq = x**2 + y**2 + z**2
        
        I[0, 0] += m * (y**2 + z**2)
        I[1, 1] += m * (x**2 + z**2)
        I[2, 2] += m * (x**2 + y**2)
        
        I[0, 1] -= m * x * y
        I[0, 2] -= m * x * z
        I[1, 2] -= m * y * z
    
    # 对称化张量
    I[1, 0] = I[0, 1]
    I[2, 0] = I[0, 2]
    I[2, 1] = I[1, 2]
    
    # 对角化得到主转动惯量
    principal_moments, eigenvectors = np.linalg.eigh(I)
    
    # 按从大到小排序
    sort_idx = np.argsort(principal_moments)[::-1]
    principal_moments = principal_moments[sort_idx]
    conversion_factor = 16.85763  # amu·Å² -> cm^-1
    principal_moments_cm = principal_moments * conversion_factor
    eigenvectors = eigenvectors[:, sort_idx]

    # Molecule.inertia is a read-only property in chem/molecule.py.
    # Store calculated principal moments on a dedicated attribute.
    molecule.inertia_principal = principal_moments
    molecule.inertia_cm = principal_moments_cm
    molecule.inertia_tensor = I
    molecule.inertia_eigenvectors = eigenvectors
    return molecule

def classify_molecule_type(principal_moments: np.ndarray) -> str:
    """
    根据主转动惯量分类分子类型
    
    分类规则:
    - 线性分子: I_a ≈ 0, I_b = I_c
    - 对称陀螺: I_a ≠ I_b = I_c (扁长) 或 I_a = I_b ≠ I_c (扁圆)
    - 不对称陀螺: I_a ≠ I_b ≠ I_c
    - 球陀螺: I_a = I_b = I_c
    """
    I_a, I_b, I_c = principal_moments
    
    # 检查是否为线性分子
    if I_a < 1e-6:
        if abs(I_b - I_c) / I_b < 0.01:
            return "线性分子 (Linear molecule)"
        else:
            return "准线性分子 (Near-linear molecule)"
    
    # 检查球陀螺
    tolerance = 0.01
    if (abs(I_a - I_b) / I_a < tolerance and 
        abs(I_b - I_c) / I_b < tolerance):
        return "球陀螺 (Spherical top)"
    
    # 检查对称陀螺
    if abs(I_a - I_b) / I_a < tolerance:
        return "扁圆对称陀螺 (Oblate symmetric top)"
    if abs(I_b - I_c) / I_b < tolerance:
        return "扁长对称陀螺 (Prolate symmetric top)"
    
    return "不对称陀螺 (Asymmetric top)"

def calculate_multipole_moments(molecule, charges) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    计算分子的偶极矩、四极矩和八极矩
    
    参数:
        coordinates: 原子坐标数组, shape (n_atoms, 3), 单位: Angstrom
        charges: 原子电荷数组, shape (n_atoms,), 单位: e (elementary charge)
    
    返回:
        dipole: 偶极矩向量, shape (3,), 单位: Debye
        quadrupole: 四极矩张量, shape (3, 3), 单位: Buckingham
        octupole: 八极矩张量, shape (3, 3, 3), 单位: e * Å^3
    """
    # 确保输入格式正确
    coordinates = np.asarray(np.array(molecule.coordinates))
    charges = np.asarray(charges)
    print(f"charges: {charges}")
    
    if coordinates.shape[0] != charges.shape[0]:
        raise ValueError("坐标和电荷数组的长度不匹配")
    
    n_atoms = len(charges)
    
    # 计算中心（对于中性分子使用几何中心）
    total_charge = np.sum(charges)
    if abs(total_charge) > 1e-10:
        center = np.sum(charges[:, np.newaxis] * coordinates, axis=0) / total_charge
    else:
        center = np.mean(coordinates, axis=0)
    
    # 将坐标平移到中心
    coords_centered = coordinates - center
    
    # 转换单位: Angstrom -> Bohr (1 Å = 1.889726125 Bohr)
    angstrom_to_bohr = 1.889726125
    coords_bohr = coords_centered * angstrom_to_bohr
    
    # 计算偶极矩
    dipole_bohr = np.sum(charges[:, np.newaxis] * coords_bohr, axis=0)
    # 转换单位: e * Bohr -> Debye (1 e * Bohr = 2.541746 Debye)
    dipole = dipole_bohr * 2.541746
    
    # 计算四极矩
    quadrupole = np.zeros((3, 3))
    for i in range(n_atoms):
        q = charges[i]
        r = coords_bohr[i]
        for j in range(3):
            for k in range(3):
                quadrupole[j, k] += 0.5 * q * (3 * r[j] * r[k] - np.sum(r**2) * (j == k))
    # 单位: e * Bohr^2, 1 e * Bohr^2 = 0.529177 Debye * Å
    quadrupole = quadrupole * 0.529177
    
    # 计算八极矩
    octupole = np.zeros((3, 3, 3))
    for i in range(n_atoms):
        q = charges[i]
        r = coords_bohr[i]
        r_sq = np.sum(r**2)
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    octupole[j, k, l] += q * (5 * r[j] * r[k] * r[l] - 
                                              r_sq * (r[j] * (k == l) + 
                                                     r[k] * (j == l) + 
                                                     r[l] * (j == k)))
    # 单位: e * Bohr^3
    octupole = octupole * (angstrom_to_bohr**3)
        # 计算四极矩的主分量
    quad_eigvals = np.linalg.eigvalsh(quadrupole)
    
    molecule.dipole = np.linalg.norm(dipole)
    molecule.dipole_moment = dipole
    molecule.quadrupole = list(quad_eigvals)
    molecule.quadrupole_moment = quadrupole
    molecule.octupole_moment = octupole

    #print(f"  |μ| = {np.linalg.norm(dipole):.6f}")
    #print("\n四极矩特征值 (单位: Buckingham):")
    #print(f"  λ1 = {quad_eigvals[0]:.6f}")
    #print(f"  λ2 = {quad_eigvals[1]:.6f}")
    #print(f"  λ3 = {quad_eigvals[2]:.6f}")
    print("dipole,", dipole)
    print(" quadrupole,", quadrupole)
    print("octupole,", octupole)
    #return dipole, quadrupole, octupole
    return molecule

def find_center(molecule):
    coor = molecule.coordinates
    mass = [atom.mass for atom in molecule.Atoms]
    x = [rr[0] for rr in coor]
    y = [rr[1] for rr in coor]
    z = [rr[2] for rr in coor]

    def calc_cog(coor):
        return [np.mean(x), np.mean(y), np.mean(z)]

    def calc_com(coor, mass):
        total_mass = sum(mass)
        x_mass = [coor[ii][0] * mass[ii] / total_mass for ii in range(len(coor))]
        y_mass = [coor[ii][1] * mass[ii] / total_mass for ii in range(len(coor))]
        z_mass = [coor[ii][2] * mass[ii] / total_mass for ii in range(len(coor))]
        return [sum(x_mass), sum(y_mass), sum(z_mass)]

    def calc_cob(coor):
        return [(min(x) + max(x)) / 2.0, (min(y) + max(y)) / 2.0, (min(z) + max(z)) / 2.0]

    def calc_coor_range(coor):
        return [[max(x), min(x)], [max(y), min(y)], [max(z), min(z)]]
    
    molecule.cog = calc_cog(coor)
    molecule.com = calc_com(coor,mass)
    molecule.cob = calc_cob(coor)
    molecule.size = calc_coor_range(coor)

    #if center_type == "cog":
    #    return calc_cog(coor)
    #elif center_type == "com":
    #    return calc_com(coor, mass)
    #elif center_type == "cob":
    #    return calc_cob(coor)
    #elif center_type == "coor_range":
    #    return calc_coor_range(coor)
    return molecule
