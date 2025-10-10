

import math
import numpy as np

# 定数
R = 6378137.0
MAX_LAT = 85.05112878
PI_RANGE = R * math.pi  # 約20,037,508.342789244

def clip_lat(lat_deg):
    return max(min(lat_deg, MAX_LAT), -MAX_LAT)

def normalize_lon(lon_deg, lon0_deg=0.0):
    d = lon_deg - lon0_deg
    d = (d + 180.0) % 360.0 - 180.0
    return lon0_deg + d

def web_mercator_xy(lon_deg, lat_deg):
    lat_c = clip_lat(lat_deg)
    lon_rad = math.radians(lon_deg)
    lat_rad = math.radians(lat_c)
    x = R * lon_rad
    y = R * math.log(math.tan(math.pi/4 + lat_rad/2))
    return np.array([x, y])

def transverse_mercator_xy(lon_deg, lat_deg):
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    B = math.cos(lat) * math.sin(lon)
    x = 0.5 * R * math.log((1 + B) / (1 - B))
    y = R * math.atan(math.tan(lat) / math.cos(lon))
    return np.array([x, y])

def ruby_transform_lonlat(lon_deg, lat_deg):
    lat_r = math.radians(lat_deg)
    lon_r = math.radians(lon_deg)
    x = math.cos(lat_r) * math.cos(lon_r)
    y = math.cos(lat_r) * math.sin(lon_r)
    z = math.sin(lat_r)
    # X軸周りに -90度回転
    x_r = x
    y_r = z
    z_r = -y
    lat_new = math.degrees(math.asin(z_r))
    lon_new = math.degrees(math.atan2(y_r, x_r))
    return lon_new, lat_new

def rotate_ccw90_xy(pt):
    x, y = pt
    return np.array([-y, x])

# テスト点
points = [
    (170, 30),
    (-170, -30),
    (120, 60),
    (0, 0),
    (45, 45),
    (-45, -45)
]

print(f"{'lon':>8} {'lat':>8} | {'横メルカトル_X':>15} {'横メルカトル_Y':>15} | {'Ruby回転_X':>15} {'Ruby回転_Y':>15} | {'差分_X':>12} {'差分_Y':>12} | {'offset_X':>8} {'offset_Y':>8}")
print("-"*120)

for lon, lat in points:
    # 横メルカトル投影
    tm = transverse_mercator_xy(lon, lat)

    # Ruby風変換で新緯度経度
    lon_new, lat_new = ruby_transform_lonlat(lon, lat)
    lon_new = normalize_lon(lon_new, 0.0)

    # Ruby風変換→Webメルカトル投影
    rb_xy = web_mercator_xy(lon_new, lat_new)
    # 90度反時計回り回転
    rb_rot = rotate_ccw90_xy(rb_xy)

    # ±1倍の R·π のオフセットをX,Y軸方向に試し最小差分を探索
    best = None
    best_norm = float('inf')
    best_offset = (0, 0)

    for ox in (-1, 0, 1):
        for oy in (-1, 0, 1):
            candidate = rb_rot + np.array([ox * PI_RANGE, oy * PI_RANGE])
            diff = candidate - tm
            norm = np.linalg.norm(diff)
            if norm < best_norm:
                best_norm = norm
                best = candidate
                best_offset = (ox, oy)

    diff_final = best - tm

    print(f"{lon:8.3f} {lat:8.3f} | {tm[0]:15.6f} {tm[1]:15.6f} | {best[0]:15.6f} {best[1]:15.6f} | {diff_final[0]:12.6f} {diff_final[1]:12.6f} | {best_offset[0]:8d} {best_offset[1]:8d}")


