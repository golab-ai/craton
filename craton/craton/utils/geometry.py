import itertools
import math
import sys

import numpy as np
from numpy.linalg import norm
from scipy.spatial.transform import Rotation

epsilon = 1e-7

import math

# 计算点到线的距离
def point_to_line_distance(point, line_start, line_end):
    x, y = point
    x1, y1 = line_start
    x2, y2 = line_end
    numerator = abs((y2-y1)*x - (x2-x1)*y + x2*y1 - y2*x1)
    denominator = math.sqrt((y2-y1)**2 + (x2-x1)**2)
    distance = numerator/denominator
    return distance

# 计算点到面的距离
def point_to_plane_distance(point, plane_normal, plane_point):
    x, y, z = point
    a, b, c = plane_normal
    x0, y0, z0 = plane_point
    numerator = abs(a*x + b*y + c*z - a*x0 - b*y0 - c*z0)
    denominator = math.sqrt(a**2 + b**2 + c**2)
    distance = numerator/denominator
    return distance

# 计算线和线的夹角
def line_to_line_angle(line1_start, line1_end, line2_start, line2_end):
    # 计算线1和线2的向量
    line1_start = np.array(line1_start)
    line1_end = np.array(line1_end)
    line2_start = np.array(line2_start)
    line2_end = np.array(line2_end)
    line1 = line1_start - line1_end
    line2 = line2_start - line2_end

    # 计算两条向量的点积
    dot_product = np.dot(line1, line2)

    # 计算两条向量的长度
    length1 = np.linalg.norm(line1)
    length2 = np.linalg.norm(line2)
    # 计算两条线之间的夹角（单位为弧度）
    #if dot_product < 0:
    #    print(dot_product) / (length1 * length2),"aaa")
    #    print(np.arccos())
    #    angle = np.pi - np.arccos(abs(dot_product) / (length1 * length2))
    #else:
    #    angle = np.arccos(dot_product / (length1 * length2))
    # 计算两条线之间的夹角（单位为弧度）
    angle = np.arccos(dot_product / (length1 * length2))
    # 将弧度转换为角度
    angle_degrees = np.degrees(angle)
    return angle_degrees

# 计算线和线的距离
def line_to_line_distance(line1_start, line1_end, line2_start, line2_end):
    x1, y1 = line1_end[0]-line1_start[0], line1_end[1]-line1_start[1]
    x2, y2 = line2_end[0]-line2_start[0], line2_end[1]-line2_start[1]
    x0, y0 = line1_start[0]-line2_start[0], line1_start[1]-line2_start[1]
    numerator = abs(x1*y2 - y1*x2 + x0*y2 - y0*x2)
    denominator = math.sqrt(x1**2 + y1**2)
    distance = numerator/denominator
    return distance

# 计算线和面的距离
def line_to_plane_distance(line_start, line_end, plane_normal, plane_point):
    x1, y1, z1 = line_start
    x2, y2, z2 = line_end
    a, b, c = plane_normal
    x0, y0, z0 = plane_point
    numerator = abs(a*x1 + b*y1 + c*z1 - a*x0 - b*y0 - c*z0)
    denominator = math.sqrt(a**2 + b**2 + c**2)
    distance = numerator/denominator
    return distance

# 计算线和面的夹角
def line_to_plane_angle(line_start, line_end, plane_normal):
    x1, y1, z1 = line_start
    x2, y2, z2 = line_end
    a, b, c = plane_normal
    x, y, z = x2-x1, y2-y1, z2-z1
    numerator = abs(a*x + b*y + c*z)
    denominator = math.sqrt(a**2 + b**2 + c**2) * math.sqrt(x**2 + y**2 + z**2)
    cos_theta = numerator / denominator
    angle = math.degrees(math.acos(cos_theta))
    return angle

def normal_vector_plane(point1, point2, point3):
    # 计算向量 v1 和 v2
    point1 = np.array(point1)
    point2 = np.array(point2)
    point3 = np.array(point3)
    v1 = point2 - point1
    v2 = point3 - point1
    #v1 = [point2[0] - point1[0], point2[1] - point1[1], point2[2] - point1[2]]
    #v2 = [point3[0] - point1[0], point3[1] - point1[1], point3[2] - point1[2]]
    
    # 计算叉积
    #cross_product = [
    #    v1[1]*v2[2] - v1[2]*v2[1],
    #    v1[2]*v2[0] - v1[0]*v2[2],
    #    v1[0]*v2[1] - v1[1]*v2[0]
    #]
    cross_product = np.cross(v1,v2)
    length = np.linalg.norm(cross_product)
    normal_vector = cross_product / length

    # 计算长度
    #length = math.sqrt(cross_product[0]**2 + cross_product[1]**2 + cross_product[2]**2)
    
    # 归一化向量
    #normal_vector = [cross_product[0]/length, cross_product[1]/length, cross_product[2]/length]
    
    return normal_vector

def plane_to_plane_angle(plane1,plane2,normal_vector=False,cross_line = False):
    n1 = normal_vector_plane(*plane1)
    n2 = normal_vector_plane(*plane2)
    dot_product = np.dot(n1,n2)
    length1 = np.linalg.norm(n1)
    length2 = np.linalg.norm(n2)
    cosine = dot_product / length1 / length2
    if cosine > 1.0:
        cosine = 1.0
    elif cosine < -1.0:
        cosine = -1.0
    angle = math.degrees(math.acos(cosine))
    if cross_line:
        cl = np.cross(n1,n2)
        if normal_vector:
            return angle, cl, n1, n2
        else:
            return angle, cl
    else:
        if normal_vector:
            return angle, n1, n2
        else:
            return angle

def distance_between_planes(plane1, plane2):
    n1, p1 = plane1
    n2, p2 = plane2
    distance = abs((p2[0] - p1[0]) * n1[0] + (p2[1] - p1[1]) * n1[1] + (p2[2] - p1[2]) * n1[2]) / math.sqrt(n1[0]**2 + n1[1]**2 + n1[2]**2)
    return distance

def rotate_about_center(points, center, angle):
    # 将所有点沿着中心点平移
    shifted_points = points - center
    
    # 将角度转换为弧度
    angle_rad = math.radians(angle)
    
    # 计算旋转矩阵
    rotation_matrix = np.array([[math.cos(angle_rad), -math.sin(angle_rad)],
                                [math.sin(angle_rad), math.cos(angle_rad)]])
    
    # 对所有点进行旋转
    rotated_points = np.dot(shifted_points, rotation_matrix)
    
    # 将所有点沿着中心点的反方向平移
    final_points = rotated_points + center
    
    return final_points

def rotate_about_axis(points, axis, angle):
    # 将所有点沿着轴平移
    shifted_points = points - axis
    
    # 将角度转换为弧度
    angle_rad = math.radians(angle)
    
    # 计算旋转矩阵
    cos_theta = math.cos(angle_rad)
    sin_theta = math.sin(angle_rad)
    ux, uy, uz = axis / np.linalg.norm(axis)
    rotation_matrix = np.array([[cos_theta + ux**2*(1-cos_theta), ux*uy*(1-cos_theta)-uz*sin_theta, ux*uz*(1-cos_theta)+uy*sin_theta],
                                [uy*ux*(1-cos_theta)+uz*sin_theta, cos_theta+uy**2*(1-cos_theta), uy*uz*(1-cos_theta)-ux*sin_theta],
                                [uz*ux*(1-cos_theta)-uy*sin_theta, uz*uy*(1-cos_theta)+ux*sin_theta, cos_theta+uz**2*(1-cos_theta)]])
    
    # 对所有点进行旋转
    rotated_points = np.dot(shifted_points, rotation_matrix)
    
    # 将所有点沿着轴的反方向平移
    final_points = rotated_points + axis
    
    return final_points

def calculate_axis_direction(point1, point2):
    # 计算轴的方向向量
    axis_direction = np.array([point2[0] - point1[0], point2[1] - point1[1], point2[2] - point1[2]])
    
    return axis_direction

def clockwise_or_counterclockwise(a, b, c, d):
    angle, cl = plane_to_plane_angle([a,b,c],[[0,1,0],[0,0,0],[1,0,0]],cross_line=True)
    angle = 180 - angle
    #if angle >= 90:
    #    angle = 180 - angle
    print(-angle)
    a = np.array(rotation_dihedral([0,0,0],cl,a,-angle))
    b = np.array(rotation_dihedral([0,0,0],cl,b,-angle))
    c = np.array(rotation_dihedral([0,0,0],cl,c,-angle))
    d = np.array(rotation_dihedral([0,0,0],cl,d,-angle))
    clock = np.cross(b[:2]-a[:2],c[:2]-b[:2])
    if d[2] >= 0:
        if clock > 0:
            cc = "顺"
        else:
            cc = "逆"
    else:
        if clock > 0:
            cc = "逆"
        else:
            cc = "顺"
    return a,b,c,d,clock,cc,angle,cl

def calculate_distance(a, b):
    
    # calculate the bond length of a and b atoms
    r = ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2) ** 0.5
    return r

def calculate_angle(a, b, c):
    # calculate the angle of a,b,c atoms
    sr1 = (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2
    sr2 = (c[0] - b[0]) ** 2 + (c[1] - b[1]) ** 2 + (c[2] - b[2]) ** 2
    sr3 = (c[0] - a[0]) ** 2 + (c[1] - a[1]) ** 2 + (c[2] - a[2]) ** 2
    cos_value = (sr1 + sr2 - sr3) / sr1**0.5 / sr2**0.5 / 2.0
    cos_value = np.clip(cos_value, -1.0, 1.0)
    tht = math.acos(cos_value)
    tht1 = tht * 180.0 / 3.1415926
    return tht1

def calculate_dihedral(a,b,c,d):
    angle, n1, n2 = plane_to_plane_angle([a,b,c],[b,c,d],normal_vector=True)
    rad = np.array(d) - np.array(a)
    if np.dot(n1,rad) < 0.0:
        return angle * -1
    return angle

def calculate_dihedral_old(a, b, c, d):
    # calculate the dihedral angle of a,b,c,d atoms
    e = [(b[0] + c[0]) / 2.0, (b[1] + c[1]) / 2.0, (b[2] + c[2]) / 2.0]
    f = [(a[0] + d[0]) / 2.0, (a[1] + d[1]) / 2.0, (a[2] + d[2]) / 2.0]

    rba = [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
    rbc = [c[0] - b[0], c[1] - b[1], c[2] - b[2]]
    rbd = [d[0] - b[0], d[1] - b[1], d[2] - b[2]]
    rad = [d[0] - a[0], d[1] - a[1], d[2] - a[2]]
    ref = [f[0] - e[0], f[1] - e[1], f[2] - e[2]]

    nabc = [rba[1] * rbc[2] - rbc[1] * rba[2], rba[2] * rbc[0] - rbc[2] * rba[0], rba[0] * rbc[1] - rbc[0] * rba[1]]
    nbcd = [rbc[1] * rbd[2] - rbd[1] * rbc[2], rbc[2] * rbd[0] - rbd[2] * rbc[0], rbc[0] * rbd[1] - rbd[0] * rbc[1]]
    dot_produce_nabcd = nabc[0] * nbcd[0] + nabc[1] * nbcd[1] + nabc[2] * nbcd[2]
    dot_produce_nabc = nabc[0] * nabc[0] + nabc[1] * nabc[1] + nabc[2] * nabc[2]
    dot_produce_nbcd = nbcd[0] * nbcd[0] + nbcd[1] * nbcd[1] + nbcd[2] * nbcd[2]
    cos_value = dot_produce_nabcd / (dot_produce_nabc * dot_produce_nbcd) ** 0.5
    if cos_value > 1.0:
        cos_value = 1.0
    elif cos_value < -1.0:
        cos_value = -1.0
    tht = math.acos(cos_value)
    tht1 = tht * 180.0 / 3.1415926

    dir1 = ref[0] * nabc[0] + ref[1] * nabc[1] + ref[2] * nabc[2]
    dir2 = ref[0] * nbcd[0] + ref[1] * nbcd[1] + ref[2] * nbcd[2]
    dir3 = nabc[0] * rad[0] + nabc[1] * rad[1] + nabc[2] * rad[2]
    if abs(dir3) < 0.00005:
        dir3 = 0.0
    if dir1 * dir2 >= 0:
        tht1 = 180.0 - tht1
    if dir3 > 0.0:
        tht1 = tht1 * -1
    return tht1

def calc_stru_para(this_coord):
    if len(this_coord) == 2:
        value = calculate_distance(this_coord[0], this_coord[1])
    elif len(this_coord) == 3:
        value = calculate_angle(this_coord[0], this_coord[1], this_coord[2])
    elif len(this_coord) == 4:
        value = calculate_dihedral(this_coord[0], this_coord[1], this_coord[2], this_coord[3])
    return value

def rotation_dihedral(b, c, d, tht, radian=False):
    """
    Rotate the atom d along the axis b-c by magnitude of angle. The rotation is
    anti-clockwise if angle > 0.

    Args:
        b: [float, float, float], cartesian coordinate of the atom in rotation
            axis
        c: [float, float, float], cartesian coordinate of the atom in rotation
            axis
        d: [float, float, float], cartesian coordinate of the atom to be rotated
        tht (float): value of angle to be changed, default unit: degrees
        radian (boolean): True when the unit of input angle is in radial

    Returns:
        [float, float, float], rotated coordiante of atom d
    """
    bc = [c[0] - b[0], c[1] - b[1], c[2] - b[2]]
    r = (bc[0] * bc[0] + bc[1] * bc[1] + bc[2] * bc[2]) ** 0.5
    nbc = [bc[0] / r, bc[1] / r, bc[2] / r]
    if not radian:
        tht = tht / 180 * np.pi
    a1 = math.cos(tht)
    a2 = 1 - a1
    a3 = math.sin(tht)
    nx = nbc[0]
    ny = nbc[1]
    nz = nbc[2]
    dx = d[0] - b[0]
    dy = d[1] - b[1]
    dz = d[2] - b[2]
    m = [
        [a1 + a2 * nx * nx, a2 * nx * ny - a3 * nz, a2 * nx * nz + a3 * ny],
        [a2 * ny * nx + a3 * nz, a1 + a2 * ny * ny, a2 * ny * nz - a3 * nx],
        [a2 * nz * nx - a3 * ny, a2 * nz * ny + a3 * nx, a1 + a2 * nz * nz],
    ]
    ddx = m[0][0] * dx + m[0][1] * dy + m[0][2] * dz + b[0]
    ddy = m[1][0] * dx + m[1][1] * dy + m[1][2] * dz + b[1]
    ddz = m[2][0] * dx + m[2][1] * dy + m[2][2] * dz + b[2]
    dd = [ddx, ddy, ddz]
    return dd

def old_change_bond(b, c, dis):
    """
    Output the translation vector on the coodinate of c to realize the bond
    length change as dis. The component of translation vector is propotional to
    the bond vector b-c. Extend the bond if dis > 0 and contract if dis < 0

    Args:
        b: [float, float, float], cartesian coordinate of the fixed atom
        c: [float, float, float], cartesian coordinate of the atom to be translated

    Returns:
        [float, float, float], translate vector on atom b
    """
    cb = [m - n for m, n in zip(c, b)]
    r = (cb[0] * cb[0] + cb[1] * cb[1] + cb[2] * cb[2]) ** 0.5
    return [x * dis / r for x in cb]

def rotate_axis_of_angle(b, c, d):
    """
    Find the rotate vector perpendicular to the plane of atom b, c and d. The
    rotate vector is calculated from cross product of angle vectors and thus
    positive rotate angles around this vector will expand the bond angle. In
    the calculation the system is translated such that the center atom c is
    placed at origin. Therefore the atoms to be rotated also need to be
    translated to origin as c and translated back after rotation.

    Args:
        b: [float, float, float], cartesian coordinate of an atom in the angle
        c: [float, float, float], cartesian coordinate of the centering atom in
            the angle
        d: [float, float, float], cartesian coordinate of an atom in the angle

    Returns:
        [float, float, float], normalized rotation axis through [0, 0, 0]
    """
    b = np.array(b)
    c = np.array(c)
    d = np.array(d)
    vector1 = b - c
    vector2 = d - c
    rotation_axis = np.cross(vector1, vector2)
    r = np.linalg.norm(rotation_axis)
    if r < epsilon:
        raise ValueError("Rotate input angle equals 180 degrees with tolerance %s" % epsilon)
    return rotation_axis / r

def old_change_angle(b, c, d, angle, radian=False):
    """
    Change the bond angle of atom b, c and d by magnitude of angle. The atoms b
    and c will be fixed while atom d is rotated. The rotate vector is
    calculated from cross product of angle vectors and thus positive rotate
    angles around this vector will expand the bond angle.

    Args:
        b: [float, float, float], cartesian coordinate of an atom in the angle
        c: [float, float, float], cartesian coordinate of the centering atom in
            the angle
        d: [float, float, float], cartesian coordinate of an atom in the angle
        angle (float): value of angle to be changed, default unit: degrees
        radian (boolean): True when the unit of input angle is in radial

    Returns:
        [float, float, float], rotated coordiante of atom d
    """
    b = np.array(b)
    c = np.array(c)
    d = np.array(d)

    if not radian:
        angle = angle / 180 * np.pi
    rotate_axis = rotate_axis_of_angle(b, c, d)
    rotater = Rotation.from_rotvec(rotate_axis * angle)
    d -= c
    d = rotater.apply(d)
    d += c
    return list(d)

def change_atom_coord(b, c, tht, change_atom, atom_arr, coord_dict):
    for d in change_atom:
        dd = rotation_dihedral(coord_dict[b], coord_dict[c], coord_dict[d], tht)
        atom_arr[int(d) - 1][4] = "%12.5f" % dd[0]
        atom_arr[int(d) - 1][5] = "%12.5f" % dd[1]
        atom_arr[int(d) - 1][6] = "%12.5f" % dd[2]
    return atom_arr

def calc_radius(coors):
    centroids = np.array(find_center(coors))
    max_dist = 0
    for coor in coors:
        max_dist = max(max_dist, np.linalg.norm(np.array(coor) - centroids))
    return max_dist

def find_center(coor, mass=None, center_type="cog"):
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

    if center_type == "cog":
        return calc_cog(coor)
    elif center_type == "com":
        return calc_com(coor, mass)
    elif center_type == "cob":
        return calc_cob(coor)
    elif center_type == "coor_range":
        return calc_coor_range(coor)

def  change_bond(a,b,p,d):
    d = -1 * d
    # 定义向量
    v = np.array([a[0]-b[0],a[1]-b[1],a[2]-b[2]])

    # 定义要移动的点的坐标
    p = np.array(p)

    # 定义要移动的距离

    # 计算移动的位移向量
    u = v / np.linalg.norm(v) * d

    # 将位移向量添加到点的坐标上
    p_new = p + u

    # 输出新的点的坐标
    return list(p_new)

def cross_product_plane(a,b,c):
    vec_AC = [b[0] - a[0], b[1] - a[1], b[2] - a[2]]
    vec_AB = [b[0] - c[0], b[1] - c[1], b[2] - c[2]]

    # 计算叉积
    cross_product = [vec_AC[1]*vec_AB[2] - vec_AC[2]*vec_AB[1],
                     vec_AC[2]*vec_AB[0] - vec_AC[0]*vec_AB[2],
                     vec_AC[0]*vec_AB[1] - vec_AC[1]*vec_AB[0]]
    return cross_product

def change_angle(a,b,c,p,angle):
    cross_product = cross_product_plane(a,b,c)
    
    v = [(b[ii] - cross_product[ii]) for ii in range(3)]
    return rotation_dihedral(v,b,p,angle)

def change_dihedral(a,b,c,d,p,angle):
    return rotation_dihedral(b,c,p,angle)

def rotate_to_axis(v,axis,coordinates):
    """
    将向量v旋转到与Z轴平行
    
    参数:
    v -- 输入向量(3维)
    
    返回:
    R -- 旋转矩阵(3x3)
    v_rotated -- 旋转后的向量(R @ v)
    """
    coordinates = np.array(coordinates)
    z_axis = np.array([0 if ii != axis else 1 for ii in range(3)])
    v = np.array(v, dtype=float)
    v = v / norm(v)
    
    
    rot_axis = np.cross(v, z_axis)
    if norm(rot_axis) < 1e-10:  # 已经平行于Z轴
        if v[2] > 0:  # 已经指向Z轴正方向
            return np.eye(3), v
        else:  # 指向Z轴反方向
            return -np.eye(3), -v
    
    rot_axis = rot_axis / norm(rot_axis)
    rot_angle = np.arccos(np.dot(v, z_axis))
    
    # 罗德里格斯旋转公式
    K = np.array([
        [0, -rot_axis[2], rot_axis[1]],
        [rot_axis[2], 0, -rot_axis[0]],
        [-rot_axis[1], rot_axis[0], 0]
    ])
    R = np.eye(3) + np.sin(rot_angle) * K + (1 - np.cos(rot_angle)) * (K @ K)
    
    v_rotated = R @ v
    return ( R @ coordinates.T).T # 转置后相乘再转置回来
    return R, v_rotated

def tetrahedron(center,distance,vertexs,):
    def center_distance():
        #distance = np.sqrt(sum([(center[ii] - vertex[ii])**2 for ii in range(3)]))
        A = [center[0],center[1],center[2] + distance] 
        B = [center[0] - 0.943 * distance,center[1],center[2] - 0.333*distance]
        C = [center[0] + 0.471 * distance,center[1] + 0.816 * distance,center[2] - 0.333*distance]
        D = [center[0] + 0.471 * distance,center[1] - 0.816 * distance,center[2] - 0.333*distance]
        return A,B,C,D

    def center_one_vertex():
        
        if v[0]**2+v[1]**2 == 0.0:
            A0 = [A[0], A[1] + abs(length * v[2])/np.sqrt(v[2]**2+v[1]**2),A[2] - abs(length * v[1])/np.sqrt(v[2]**2+v[1]**2)]
        else:
            if v[0] * v[1] >= 0.0:
                A0 = [A[0] + abs(length * v[1])/np.sqrt(v[0]**2+v[1]**2), A[1] - abs(length * v[0])/np.sqrt(v[0]**2+v[1]**2),A[2]]
            else:
                A0 = [A[0] + abs(length * v[1])/np.sqrt(v[0]**2+v[1]**2), A[1] + abs(length * v[0])/np.sqrt(v[0]**2+v[1]**2),A[2]]
        A1 = rotation_dihedral(center,A,A0,120)
        A2 = rotation_dihedral(center,A,A0,240)
        #A0 = rotation_dihedral(vertex,center,A0,60)
        #A1 = rotation_dihedral(vertex,center,A1,60)
        #A2 = rotation_dihedral(vertex,center,A2,60)
        return A0, A1, A2

    def center_two_vertexs():
        A1 = rotation_dihedral(center,A,A0,120)
        A2 = rotation_dihedral(center,A,A0,240)
        return A1, A2

    def center_three_vertexs():
        A20 = rotation_dihedral(center,A,A0,120)
        A21 = rotation_dihedral(center,A,A0,240)
        angle1 = line_to_line_angle(A1,A,A20,A)
        angle2 = line_to_line_angle(A1,A,A21,A)
        #return [A20 if angle1 < angle2 else A21]
        return [A20 if angle1 > angle2 else A21]

    if center is not None:
        if len(vertexs) == 0:
            if distance is None:
                return 
            else:
                return center_distance()
        else:
            height = distance * np.sin(np.deg2rad(19.4712))
            length = distance * np.cos(np.deg2rad(19.4712))
            #A = change_bond(vertexs[0],center,center,-height)
            A = change_bond(vertexs[0],center,center,height)
            v = [(center[0]-A[0])/height,(center[1]-A[1])/height,(center[2]-A[2])/height]
            if len(vertexs) == 1:
                return center_one_vertex()
            elif len(vertexs) == 2:
                A0 = vertexs[1]
                return center_two_vertexs()
            elif len(vertexs) == 3:
                A0 = vertexs[1]
                A1 = vertexs[2]
                return center_three_vertexs()
    else:
        if len(vertexs) == 2:
            pass

def triangle(center,distance,vertexs,):
    def center_distance():
        vertex1 = [center[0] - distance, center[1],center[2]]
        vertex2 = [center[0] + 0.5*distance, center[1] + np.sqrt(3/4)*distance ,center[2]]
        vertex3 = [center[0] + 0.5*distance, center[1] - np.sqrt(3/4)*distance ,center[2]]
        return vertex1, vertex2, vertex3

    def center_one_vertex():
        for ii in range(3):
            if center[ii] - vertex1[ii] != 0.0:
                cc = ii
                break
        z = sum([(center[ii] - vertex1[ii])*(center[ii] - 1)/(center[cc]-vertex1[cc]) if ii != cc 
                 else  (center[ii] - vertex1[ii])*(center[ii])/(center[cc]-vertex1[cc]) for ii in range(3)])
        v = [z if ii == cc else 1 for ii in range(3)]
        vertex2 = rotation_dihedral(v,center,vertex1,120)
        vertex3 = rotation_dihedral(v,center,vertex1,-120)
        return vertex2, vertex3

    def center_two_vertexs():
        vertex3 = rotation_dihedral(vertex1,center,vertex2,180)
        return [vertex3]
    
    if center is not None:
        if len(vertexs) == 0:
            if distance is None:
                return 
            else:
                return center_distance()
        else:
            if len(vertexs) == 1:
                vertex1 = vertexs[0]
                return center_one_vertex()
            elif len(vertexs) == 2:
                vertex1 = vertexs[0]
                vertex2 = vertexs[1]
                return center_two_vertexs()
    else:
        if len(vertexs) == 2:
            pass

def line(center,distance,vertexs,):
    def center_distance():
        vertex1 = [center[0],center[1],center[2] + distance]
        vertex2 = [center[0],center[1],center[2] - distance]
        return vertex1, vertex2
    def center_one_vertex():
        vertex2 = change_bond(vertexs[0],center,center,-distance)
        return [vertex2]
    
    if center is not None:
        if len(vertexs) == 0:
            if distance is None:
                return 
            else:
                return center_distance()
        else:
            return center_one_vertex()

def center_of_ring(points):
    a = np.array(points[0])
    b = np.array(points[1])
    c = np.array(points[2])
    ra2 = np.linalg.norm(a)
    rb2 = np.linalg.norm(b)
    rc2 = np.linalg.norm(c)
    Vba = b - a
    Vca = c - a
    D2 = ra2 - rb2
    D3 = ra2 - rc2
    A1 = a[1]*b[2] - a[1]*c[2] - a[2]*b[1] + a[2]*c[1] + b[1]*c[2] - c[1]*b[2]
    B1 = -a[0]*b[2] + a[0]*c[2] + a[2]*b[0] - a[2]*c[0] - b[0]*c[2] + c[0]*b[2]
    C1 = a[0]*b[1] - a[0]*c[1] - a[1]*b[0] + a[1]*c[0] + b[0]*c[1] - c[0]*b[1]
    D1 = -a[0]*b[1]*c[2] + a[0]*c[1]*b[2] + b[0]*a[1]*c[2] - c[0]*a[1]*b[2] - b[0]*c[1]*a[2] + c[0]*b[1]*a[2]
    M = np.array([np.array([A1,B1,C1]),Vba,Vca])
    M = np.matrix(M)
    N = np.array([D1,D2,D3])
    N = np.matrix(N)
    center = np.matmul(-M.I,N)
    return center

def square_2d(n,center,vertexs,distance,other):
    def _equ_distances(dists):
        for i in range(1,len(dists)):
            if abs(dists[0] - dists[i]) > 0.01:
                return False
        return True

    if not isinstance(n,int):
        return None
    angle = 360.0/n

    if center is None:
        if vertexs is None:
            center = [0.0,0.0,0.0]
            if distance is None:
                distance = 1.0
            vertexs =[[distance,0.0,0.0]]
        else:
            if len(vertexs) == 1:
                if distance is None:
                    distance = 1.0
                if other is None:
                    center = [vertexs[0][0] + distance,vertexs[0][1],vertexs[0][2]]
                else:
                    center = change_bond(other,vertexs[0],vertexs[0],distance)
            elif len(vertexs) == 2:
                distance_tmp = (calculate_distance(*vertexs)/(2*(1-math.cos(math.radians(angle)))))**0.5
                if distance is not None:
                    if abs(distance_tmp - distance) > 0.01:
                        return "distance error"
                distance = distance_tmp
                
            else:
                center = center_of_ring(vertexs[:3])
                distance_tmps = [calculate_distance(center,point) for point in vertexs]
                if not _equ_distances(distance_tmps):
                    return "distance error"
                if distance is not None:
                    if abs(distance_tmps[0] - distance) > 0.01:
                        return "distance error"
                distance = distance_tmps[0]
    if vertexs is None:
        if distance is None:
            distance = 1.0
        vertexs = [[center[0]+distance,center[1],center[2]]]
    


    def center_distance():
        vertex1 = [center[0] - distance, center[1],center[2]]
        vertex2 = [center[0] + distance, center[1],center[2]]
        vertex3 = [center[0], center[1] - distance,center[2]]
        vertex4 = [center[0], center[1] + distance,center[2]]
        return vertex1, vertex2, vertex3, vertex4
    
    def center_one_vertex():
        cn = -1
        for ii in range(3):
            if round(center[ii],3) - round(vertexs[0][ii],3) == 0.0:
                cn = ii
        if cn != -1:
            v = [center[ii] if ii != cn else center[ii] + 1.0 for ii in range(3)]
        else:
            cc = [ii for ii in range(3) if ii != cn][0]
            z = sum([(center[ii] - vertexs[0][ii])*(center[ii] - 1)/(center[cc]-vertexs[0][cc]) if ii != cc 
                 else  (center[ii] - vertexs[0][ii])*(center[ii])/(center[cc]-vertexs[0][cc]) for ii in range(3)])
            v = [z if ii == cc else 1 for ii in range(3)]

        vertex2 = rotation_dihedral(v,center,vertexs[0],90)
        vertex3 = rotation_dihedral(v,center,vertexs[0],180)
        vertex4 = rotation_dihedral(v,center,vertexs[0],270)

    def center_multi_vertexs():
        v = normal_vector_plane(center,vertexs[0],vertexs[1])
        v = np.array(center) - v
        exist_dihes = [calculate_dihedral(vertexs[0],center,v,vertex) for vertex  in vertexs[1:]]
        exist_dihes = [int((dihe+4)/angle) for dihe in exist_dihes]
        create_edge = [i for i in range(1,n) if i not in exist_dihes]
        create_vertexs = [rotation_dihedral(v,center,vertexs[0],ii*angle) for ii in create_edge]
        return create_vertexs