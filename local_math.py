# -*- coding: utf-8 -*-
import numpy as np
from Bio.PDB import Vector
def calculate_distance(coord1, coord2):
    return float( np.linalg.norm(coord2 - coord1))

def calculate_angle(coord1, coord2, coord3):
    vector1 = coord1 - coord2
    vector2 = coord3 - coord2
    cosine_angle = np.dot(vector1, vector2) / (np.linalg.norm(vector1) * np.linalg.norm(vector2))
    angle = np.arccos(np.clip(cosine_angle, -1.0, 1.0))
    return float( np.degrees(angle))

def calculate_dihedral(coord1, coord2, coord3, coord4):
    vector1 = coord2 - coord1
    vector2 = coord3 - coord2
    vector3 = coord4 - coord3

    normal1 = np.cross(vector1, vector2)
    normal2 = np.cross(vector2, vector3)

    cosine_angle = np.dot(normal1, normal2) / (np.linalg.norm(normal1) * np.linalg.norm(normal2))
    angle = np.arccos(np.clip(cosine_angle, -1.0, 1.0))
    # 确定二面角的符号
    if np.dot(vector2, np.cross(normal1, normal2)) < 0:
        angle = -angle

    return float(np.degrees(angle))

def angle_between_vectors(a: np.ndarray, b: np.ndarray) -> float:
    """返回两个向量之间的夹角（单位：度）"""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)

    dot = np.dot(a, b)
    norm_a = np.linalg.norm(a)
    norm_b = np.linalg.norm(b)

    if norm_a == 0 or norm_b == 0:
        raise ValueError("不能计算零向量的夹角")

    cos_theta = np.clip(dot / (norm_a * norm_b), -1.0, 1.0)  # 防止浮点误差
    angle_rad = np.arccos(cos_theta)
    angle_deg = np.degrees(angle_rad)
    return angle_deg

def dihedral_among_vectors(v1, v2, v3):
    """
    计算三个首尾相连向量 v1→v2→v3 的二面角（单位：度）
    返回 [-180, 180] 的有符号角，符号表示旋转方向
    """
    v1, v2, v3 = map(np.asarray, (v1, v2, v3))

    # 平面法向量
    n1 = np.cross(v1, v2)
    n2 = np.cross(v2, v3)

    # 归一化
    n1 = n1 / (np.linalg.norm(n1) + 1e-12)
    n2 = n2 / (np.linalg.norm(n2) + 1e-12)

    # 计算夹角
    cos_phi = np.clip(np.dot(n1, n2), -1.0, 1.0)
    angle = np.degrees(np.arccos(cos_phi))

    # 判断方向：v2 方向看，n1 到 n2 是顺时针还是逆时针
    sign = np.sign(np.dot(np.cross(n1, n2), v2))
    return angle * sign

def mix_z_matrix_to_cartesian(r1, r2, r3, r41, theta412, phi4123):
    """
    Convert Z-matrix coordinates to Cartesian coordinates.

    Parameters:
    r1, r2, r3: Cartesian coordinates of the reference atoms (numpy arrays of shape (3,))
    r41: Distance from atom 4 to atom 1
    theta412: Angle between atoms 4, 1, and 2 (in degrees)
    phi4123: Dihedral angle between atoms 4, 1, 2, and 3 (in degrees)

    Returns:
    r4: Cartesian coordinates of atom 4 (numpy array of shape (3,))
    """
    # Convert angles from degrees to radians
    theta412_rad = np.radians(theta412)
    phi4123_rad = np.radians(phi4123)

    # Convert to numpy arrays
    r1 = np.array(r1)
    r2 = np.array(r2)
    r3 = np.array(r3)

    # Calculate vectors
    v12 = r2 - r1  # Vector from atom 1 to atom 2
    v23 = r3 - r2  # Vector from atom 2 to atom 3

    # Normalize vectors
    u12 = v12 / np.linalg.norm(v12)
    u23 = v23 / np.linalg.norm(v23)

    # Calculate normal vector (perpendicular to plane formed by atoms 1, 2, and 3)
    n = np.cross(u12, u23)
    n = n / np.linalg.norm(n)

    # Calculate the vector perpendicular to u12 and n (in the plane of atoms 1, 2, and 3)
    u_perp = np.cross(n, u12)
    u_perp = u_perp / np.linalg.norm(u_perp)

    # Calculate the coordinates of atom 4
    r4 = r1 + r41 * (
            u12 * np.cos(theta412_rad) +
            u_perp * np.sin(theta412_rad) * np.cos(phi4123_rad) +
            n * np.sin(theta412_rad) * np.sin(phi4123_rad)
    )

    return r4


def calculate_third_point(x_bond, y_bond, x_angle, y_angle, bond, angle):
    """
    计算第三点的坐标。

    参数:
    - x1, y1: 第一个点的坐标
    - x2, y2: 第二个点的坐标
    - phi: 角度（单位为度）
    - d23: 点3与点2的距离

    返回:
    - (x3, y3): 第三个点的坐标
    """
    # 将角度转换为弧度
    phi_rad = np.radians(angle)

    # 计算向量 v12
    v12 = np.array([x_angle - x_bond, y_angle - y_bond])

    # 计算向量 v12 的长度
    d12 = np.linalg.norm(v12)

    # 计算单位向量 u12
    u12 = v12 / d12

    # 计算垂直于 u12 的单位向量 u12_perp
    u12_perp = np.array([-u12[1], u12[0]])

    # 计算向量 v23
    v23 = bond * (np.cos(phi_rad) * u12 + np.sin(phi_rad) * u12_perp)

    # 计算点3的坐标
    x3 = x_angle + v23[0]
    y3 = y_angle + v23[1]

    return x3, y3

def read_pdb_as_cartesian(pdb_file):
    coordinates = []
    with open(pdb_file, 'r') as file:
        for line in file:
            atom = []
            if line.startswith('ATOM'):
                atom_name = line[12:16].strip()
                residue_serial = line[22:26].strip()
                residue_name = line[17:20].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                atom.append(atom_name)
                atom.append(residue_name)
                atom.append(residue_serial)
                atom.append(x)
                atom.append(y)
                atom.append(z)
            if len(atom) > 0:
                coordinates.append(atom)
    return coordinates

def get_normal_vector(v1,v2,v3):
    return np.cross(v1-v2,v1-v3)

def mirror_about_plane(structure, plane_normal: Vector, plane_point: Vector):
    """
    将 structure 中所有原子关于“过 plane_point 且以 plane_normal 为法向的平面”做镜像。
    plane_normal 无需单位化，函数内部会归一化。
    """
    if not isinstance(plane_normal, Vector):
        plane_normal = Vector(*plane_normal)
    if not isinstance(plane_point, Vector):
        plane_point = Vector(*plane_point)

    n = np.array(plane_normal.get_array(), dtype=float)
    n = n / (np.linalg.norm(n) + 1e-16)          # 单位法向
    p0 = np.array(plane_point.get_array(), dtype=float)

    # 镜像矩阵 I - 2 n⊗n
    T = np.eye(3) - 2 * np.outer(n, n)

    for res in structure:
        for atom in res.get_atoms():
            coord = np.array(atom.get_coord(), dtype=float)
            mirrored = p0 + T @ (coord - p0)
            atom.set_coord(mirrored)