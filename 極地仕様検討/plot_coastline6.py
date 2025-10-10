import math
import matplotlib.pyplot as plt
import shapefile  # pyshpライブラリ、事前にpip install pyshpが必要

# 地球の半径[m]
R = 6378137.0

# WebメルカトルのX,Yの最大範囲（地球円周の半分）
MAX_RANGE = R * math.pi  # 約20,037,508.34m

# -----------------------------------------------
# 1: 緯度経度（度）→ Webメルカトル投影(m)変換関数
def latlon_to_webmercator(lat_deg, lon_deg):
    """
    Webメルカトル投影法:
    - 緯度の範囲は±85.05112878度でクリップ（投影法の定義範囲）
    - 経度は -180〜180度
    - 出力はメートル単位のX,Y座標
    """
    lat_deg = max(min(lat_deg, 85.05112878), -85.05112878)  # 緯度範囲制限
    lat_rad = math.radians(lat_deg)  # 緯度をラジアンに変換
    lon_rad = math.radians(lon_deg)  # 経度をラジアンに変換
    x = R * lon_rad  # X座標（経度に比例）
    y = R * math.log(math.tan(math.pi/4 + lat_rad/2))  # Y座標（緯度に非線形関数）
    return x, y

# -----------------------------------------------
# 2: 緯度経度（度）→ 横メルカトル（±180度縦範囲拡張版）投影(m)
def latlon_to_transverse_mercator_extended(lat_deg, lon_deg, lon0_deg=0):
    """
    横メルカトル投影の拡張版:
    - Y座標方向を±180度の範囲で展開（通常は±90度）
    - lon0_degは中心経線（デフォルト0度）
    - 出力はメートル単位のX,Y座標
    """
    lat_rad = math.radians(lat_deg)
    lon_rad = math.radians(lon_deg)
    lon0_rad = math.radians(lon0_deg)

    # 投影で使う角度を計算（緯度方向の変形）
    angle_deg = math.degrees(math.atan2(math.tan(lat_rad), math.cos(lon_rad - lon0_rad)))

    # Bの計算（横メルカトルのX座標成分）
    B = math.cos(lat_rad) * math.sin(lon_rad - lon0_rad)

    # X座標計算（対数関数で横方向座標を得る）
    x = 0.5 * R * math.log((1 + B) / (1 - B))

    # Y座標は緯度方向の角度をラジアンに戻して距離換算
    y = R * math.radians(angle_deg)

    return x, y

# -----------------------------------------------
# 3: Rubyのtransform_latlng関数のPython版（球面座標の回転変換）
def transform_latlng(lon_deg, lat_deg):
    """
    入力: 経度lon_deg, 緯度lat_deg (度)
    処理: 3D球面座標に変換し、X軸周りに-90度回転させてから緯度経度に戻す
    出力: 回転変換後の経度, 緯度 (度)
    """
    lat_rad = math.radians(lat_deg)
    lon_rad = math.radians(lon_deg)

    # 球面直交座標に変換
    x = math.cos(lat_rad) * math.cos(lon_rad)
    y = math.cos(lat_rad) * math.sin(lon_rad)
    z = math.sin(lat_rad)

    # X軸周りに-90度回転 (x'=x, y'=z, z'=-y)
    x_r = x
    y_r = z
    z_r = -y

    # 回転後座標から緯度経度に戻す
    lat_new = math.degrees(math.asin(z_r))
    lon_new = math.degrees(math.atan2(y_r, x_r))

    return lon_new, lat_new

# -----------------------------------------------
# 線分描画時の不連続区間を分割し描画する関数
def plot_lines_filtered(lines_proj, max_x_diff=1_000_000, max_segment_length=2_000_000, color="green"):
    """
    線の点群リストを受け取り、隣接点間の
    - X座標差分の閾値(max_x_diff)
    - 線分長さの閾値(max_segment_length)
    を超える箇所で線を分割して描画。

    color指定で線色を変えられます。
    """
    for proj_pts in lines_proj:
        segments = []
        current_segment = [proj_pts[0]]

        for i in range(1, len(proj_pts)):
            x1, y1 = proj_pts[i-1]
            x2, y2 = proj_pts[i]
            dist = math.hypot(x2 - x1, y2 - y1)  # 線分長さ

            # 距離かX方向差が大きい場合は線を切る
            if abs(x2 - x1) > max_x_diff or dist > max_segment_length:
                if len(current_segment) >= 2:
                    segments.append(current_segment)
                current_segment = [proj_pts[i]]
            else:
                current_segment.append(proj_pts[i])

        if len(current_segment) >= 2:
            segments.append(current_segment)

        # 分割したセグメントごとに描画
        for seg in segments:
            xs, ys = zip(*seg)
            plt.plot(xs, ys, color=color, linewidth=0.5)

# -----------------------------------------------
# shapefile（Natural Earth coastline）読み込み
sf = shapefile.Reader("ne_10m_coastline/ne_10m_coastline.shp")

# 各投影結果を格納するリストを用意
latlon_lines = []
webmerc_lines = []
tmerc_lines = []
ruby_transform_lines = []

# shapefile中の全てのshape（ポリライン）を処理
for shape in sf.shapes():
    lons = [pt[0] for pt in shape.points]  # 経度一覧
    lats = [pt[1] for pt in shape.points]  # 緯度一覧
    latlon_pts = list(zip(lons, lats))
    latlon_lines.append(latlon_pts)

    # 2番目グラフ用：Webメルカトル投影
    webmerc_lines.append([latlon_to_webmercator(lat, lon) for lon, lat in latlon_pts])

    # 3番目グラフ用：横メルカトル拡張投影
    tmerc_lines.append([latlon_to_transverse_mercator_extended(lat, lon) for lon, lat in latlon_pts])

    # 4番目グラフ用：Ruby風変換（緯度経度）
    ruby_transform_lines.append([transform_latlng(lon, lat) for lon, lat in latlon_pts])

# -----------------------------------------------
# 1: 元の緯度経度（度）を黒で描画
plt.figure("Original Lat/Lon (degrees)")
for pts in latlon_lines:
    lons, lats = zip(*pts)
    plt.plot(lons, lats, color="black", linewidth=0.5)
plt.xlim(-180, 180)
plt.ylim(-90, 90)
plt.xlabel("Longitude (deg)")
plt.ylabel("Latitude (deg)")
plt.title("Original Geographic Coordinates")
plt.grid(True)

# -----------------------------------------------
# 2: Webメルカトル投影(m)、青線
plt.figure("Web Mercator (EPSG:3857)")
for line in webmerc_lines:
    xs, ys = zip(*line)
    plt.plot(xs, ys, color="blue", linewidth=0.5)
plt.xlim(-MAX_RANGE, MAX_RANGE)
plt.ylim(-MAX_RANGE, MAX_RANGE)
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.title("Web Mercator Projection")
plt.grid(True)
plt.gca().set_aspect('equal', adjustable='box')

# -----------------------------------------------
# 3: 横メルカトル拡張（±180°Y方向）、緑線
plt.figure("Transverse Mercator Extended ±180° Y")
plot_lines_filtered(tmerc_lines, max_x_diff=1_000_000, max_segment_length=2_000_000, color="green")
plt.xlim(-MAX_RANGE, MAX_RANGE)
plt.ylim(-MAX_RANGE, MAX_RANGE)
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.title("Transverse Mercator Projection (Extended ±180° Y)")
plt.grid(True)
plt.gca().set_aspect('equal', adjustable='box')

# -----------------------------------------------
# 4: Ruby風変換の緯度経度（度）、黒線に変更
plt.figure("Ruby-style transform_latlng")
plot_lines_filtered(ruby_transform_lines, max_x_diff=20, max_segment_length=20, color="black")
plt.xlim(-180, 180)
plt.ylim(-90, 90)
plt.xlabel("Longitude (deg)")
plt.ylabel("Latitude (deg)")
plt.title("Ruby-style transform_latlng Projection (No Jump Lines)")
plt.grid(True)

# -----------------------------------------------
# 5: 4番目の結果をWebメルカトルに投影（紫線）
ruby_webmerc_lines = []
for line in ruby_transform_lines:
    # 緯度経度→Webメルカトル(m)に変換
    ruby_webmerc_lines.append([latlon_to_webmercator(lat, lon) for lon, lat in line])

plt.figure("Ruby-style transform_latlng → Web Mercator")
plot_lines_filtered(ruby_webmerc_lines, max_x_diff=1_000_000, max_segment_length=2_000_000, color="purple")
plt.xlim(-MAX_RANGE, MAX_RANGE)
plt.ylim(-MAX_RANGE, MAX_RANGE)
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.title("Ruby-style transform_latlng Projected to Web Mercator")
plt.grid(True)
plt.gca().set_aspect('equal', adjustable='box')

# -----------------------------------------------
plt.show()
