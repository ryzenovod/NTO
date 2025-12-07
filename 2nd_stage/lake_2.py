import math

def vvod_dannyh():
    with open('input.txt', 'r') as f:
        nachalnaya_tochka = tuple(map(float, f.readline().split()))
        konechnaya_tochka = tuple(map(float, f.readline().split()))
        skorost = float(f.readline())
        napravlenie_kamery = tuple(map(float, f.readline().split()))
        kolichestvo_lodok, kolichestvo_gor = map(int, f.readline().split())
        
        lodki = []
        for _ in range(kolichestvo_lodok):
            lodki.append(tuple(map(float, f.readline().split())))
        
        gori = []
        for _ in range(kolichestvo_gor):
            gori.append(tuple(map(float, f.readline().split())))
    
    return nachalnaya_tochka, konechnaya_tochka, skorost, napravlenie_kamery, lodki, gori

def postroit_bazis_kamery(napravlenie_kamery):
    dir_x, dir_y, dir_z = napravlenie_kamery
    dlina_napravleniya = math.sqrt(dir_x * dir_x + dir_y * dir_y + dir_z * dir_z)
    
    if dlina_napravleniya < 1e-9:
        return None
    
    vpered_x = dir_x / dlina_napravleniya
    vpered_y = dir_y / dlina_napravleniya
    vpered_z = dir_z / dlina_napravleniya
    
    if abs(dir_x) + abs(dir_y) > 1e-9:
        vpravo_temp_x = -dir_y
        vpravo_temp_y = dir_x
        vpravo_temp_z = 0.0
    else:
        vpravo_temp_x = 1.0
        vpravo_temp_y = 0.0
        vpravo_temp_z = 0.0
    
    dlina_vpravo_temp = math.sqrt(vpravo_temp_x * vpravo_temp_x + 
                                  vpravo_temp_y * vpravo_temp_y + 
                                  vpravo_temp_z * vpravo_temp_z)
    vpravo_x = vpravo_temp_x / dlina_vpravo_temp
    vpravo_y = vpravo_temp_y / dlina_vpravo_temp
    vpravo_z = vpravo_temp_z / dlina_vpravo_temp
    
    vverh_x = vpered_y * vpravo_z - vpered_z * vpravo_y
    vverh_y = vpered_z * vpravo_x - vpered_x * vpravo_z
    vverh_z = vpered_x * vpravo_y - vpered_y * vpravo_x
    
    return (vpered_x, vpered_y, vpered_z), (vpravo_x, vpravo_y, vpravo_z), (vverh_x, vverh_y, vverh_z)

def predobrabotka_gor(gori):
    obrabotannye = []
    for gora_x, gora_y, vysota_gori, radius_gori in gori:
        if vysota_gori <= 1e-9 or radius_gori <= 1e-9:
            obrabotannye.append(None)
            continue
        kvadrat_otnosheniya_radiusa = (radius_gori * radius_gori) / (vysota_gori * vysota_gori)
        kvadrat_radiusa = radius_gori * radius_gori
        obrabotannye.append((gora_x, gora_y, vysota_gori, radius_gori, 
                         kvadrat_otnosheniya_radiusa, kvadrat_radiusa))
    return obrabotannye

def lodka_v_kadre(poziciya_kamery, poziciya_lodki, bazis_kamery):
    kamera_x, kamera_y, kamera_z = poziciya_kamery
    lodka_x, lodka_y, vysota_lodki = poziciya_lodki
    vpered, vpravo, vverh = bazis_kamery
    
    vpered_x, vpered_y, vpered_z = vpered
    vpravo_x, vpravo_y, vpravo_z = vpravo
    vverh_x, vverh_y, vverh_z = vverh
    
    vektor_k_lodke_x = lodka_x - kamera_x
    vektor_k_lodke_y = lodka_y - kamera_y
    vektor_k_lodke_z = vysota_lodki - kamera_z
    
    glubina = vektor_k_lodke_x * vpered_x + vektor_k_lodke_y * vpered_y + vektor_k_lodke_z * vpered_z
    if glubina <= 1e-9:
        return False
    
    gorizontalnoe_smeshchenie = vektor_k_lodke_x * vpravo_x + vektor_k_lodke_y * vpravo_y + vektor_k_lodke_z * vpravo_z
    vertikalnoe_smeshchenie = vektor_k_lodke_x * vverh_x + vektor_k_lodke_y * vverh_y + vektor_k_lodke_z * vverh_z
    
    if abs(gorizontalnoe_smeshchenie) > glubina + 1e-12:
        return False
    if abs(vertikalnoe_smeshchenie) > glubina + 1e-12:
        return False
    return True

def nayti_blizhayshuyu_tochku_na_luche_xy(nachalo_lucha_x, nachalo_lucha_y, napravlenie_lucha_x, napravlenie_lucha_y, 
                                  cel_x, cel_y):
    vektor_k_celi_x = cel_x - nachalo_lucha_x
    vektor_k_celi_y = cel_y - nachalo_lucha_y
    kvadrat_dliny_napravleniya = napravlenie_lucha_x * napravlenie_lucha_x + napravlenie_lucha_y * napravlenie_lucha_y
    
    if kvadrat_dliny_napravleniya < 1e-9:
        return nachalo_lucha_x, nachalo_lucha_y
    
    parametr_proekcii = (napravlenie_lucha_x * vektor_k_celi_x + napravlenie_lucha_y * vektor_k_celi_y) / kvadrat_dliny_napravleniya
    parametr_proekcii = max(0.0, min(1.0, parametr_proekcii))
    
    blizhayshaya_x = nachalo_lucha_x + napravlenie_lucha_x * parametr_proekcii
    blizhayshaya_y = nachalo_lucha_y + napravlenie_lucha_y * parametr_proekcii
    return blizhayshaya_x, blizhayshaya_y

def nayti_interval_po_vysote(kamera_z, napravlenie_lucha_z, vysota_gori):
    if abs(napravlenie_lucha_z) < 1e-9:
        if kamera_z < -1e-9 or kamera_z > vysota_gori + 1e-9:
            return None
        return 0.0, 1.0
    
    t_min = -1e30
    t_max = 1e30
    
    t_na_nule = -kamera_z / napravlenie_lucha_z
    if napravlenie_lucha_z > 0:
        t_min = max(t_min, t_na_nule)
    else:
        t_max = min(t_max, t_na_nule)
    
    t_na_vershine = (vysota_gori - kamera_z) / napravlenie_lucha_z
    if napravlenie_lucha_z > 0:
        t_max = min(t_max, t_na_vershine)
    else:
        t_min = max(t_min, t_na_vershine)
    
    t_z_min = max(0.0, t_min)
    t_z_max = min(1.0, t_max)
    
    if t_z_min > t_z_max + 1e-9:
        return None
    
    return t_z_min, t_z_max

def vychislit_kvadratichnuyu(t, koef_a, koef_b, koef_c):
    return (koef_a * t + koef_b) * t + koef_c

def gora_zakryvaet_lodku(poziciya_kamery, poziciya_lodki, obrabotannaya_gora):
    if obrabotannaya_gora is None:
        return False
    
    kamera_x, kamera_y, kamera_z = poziciya_kamery
    lodka_x, lodka_y, vysota_lodki = poziciya_lodki
    gora_x, gora_y, vysota_gori, radius_gori, kvadrat_otnosheniya_radiusa, kvadrat_radiusa = obrabotannaya_gora
    
    napravlenie_lucha_x = lodka_x - kamera_x
    napravlenie_lucha_y = lodka_y - kamera_y
    napravlenie_lucha_z = vysota_lodki - kamera_z
    
    if abs(napravlenie_lucha_x) < 1e-9 and abs(napravlenie_lucha_y) < 1e-9:
        kvadrat_rasstoyaniya = (kamera_x - gora_x) * (kamera_x - gora_x) + \
                      (kamera_y - gora_y) * (kamera_y - gora_y)
        if kvadrat_rasstoyaniya > kvadrat_radiusa + 1e-8:
            return False
    else:
        blizhayshaya_x, blizhayshaya_y = nayti_blizhayshuyu_tochku_na_luche_xy(
            kamera_x, kamera_y, napravlenie_lucha_x, napravlenie_lucha_y, gora_x, gora_y)
        dx_gora = blizhayshaya_x - gora_x
        dy_gora = blizhayshaya_y - gora_y
        kvadrat_rasstoyaniya = dx_gora * dx_gora + dy_gora * dy_gora
        if kvadrat_rasstoyaniya > kvadrat_radiusa + 1e-8:
            return False
    
    interval_vysoty = nayti_interval_po_vysote(kamera_z, napravlenie_lucha_z, vysota_gori)
    if interval_vysoty is None:
        return False
    
    t_z_min, t_z_max = interval_vysoty
    t_nachalo = max(t_z_min, 1e-7)
    t_konec = t_z_max
    
    if t_nachalo > t_konec + 1e-9:
        return False
    
    smeshchenie_x = kamera_x - gora_x
    smeshchenie_y = kamera_y - gora_y
    kvadrat_rasstoyaniya_kamery = smeshchenie_x * smeshchenie_x + smeshchenie_y * smeshchenie_y
    raznica_vysoty = vysota_gori - kamera_z
    
    koef_a_xy = napravlenie_lucha_x * napravlenie_lucha_x + napravlenie_lucha_y * napravlenie_lucha_y
    koef_b_xy = 2.0 * (smeshchenie_x * napravlenie_lucha_x + smeshchenie_y * napravlenie_lucha_y)
    koef_c_xy = kvadrat_rasstoyaniya_kamery
    
    koef_a_konus = kvadrat_otnosheniya_radiusa * napravlenie_lucha_z * napravlenie_lucha_z
    koef_b_konus = -2.0 * kvadrat_otnosheniya_radiusa * raznica_vysoty * napravlenie_lucha_z
    koef_c_konus = kvadrat_otnosheniya_radiusa * raznica_vysoty * raznica_vysoty
    
    koef_a = koef_a_xy - koef_a_konus
    koef_b = koef_b_xy - koef_b_konus
    koef_c = koef_c_xy - koef_c_konus
    
    if abs(koef_a) < 1e-12:
        znachenie_na_nachale = vychislit_kvadratichnuyu(t_nachalo, koef_a, koef_b, koef_c)
        znachenie_na_konce = vychislit_kvadratichnuyu(t_konec, koef_a, koef_b, koef_c)
        min_znachenie = min(znachenie_na_nachale, znachenie_na_konce)
        return min_znachenie <= 1e-9
    
    parametr_vershiny = -koef_b / (2.0 * koef_a)
    kandidaty = [t_nachalo, t_konec]
    if t_nachalo - 1e-12 <= parametr_vershiny <= t_konec + 1e-12:
        kandidaty.append(parametr_vershiny)
    
    min_znachenie = float('inf')
    for kandidat in kandidaty:
        ogranichennyy = max(t_nachalo, min(t_konec, kandidat))
        znachenie = vychislit_kvadratichnuyu(ogranichennyy, koef_a, koef_b, koef_c)
        if znachenie < min_znachenie:
            min_znachenie = znachenie
    
    return min_znachenie <= 1e-9

def reshit_staticheskiy_sluchay(poziciya_kamery, lodki, obrabotannye_gori, bazis_kamery):
    vidimye_lodki = []
    for nomer_lodki, lodka in enumerate(lodki, start=1):
        if not lodka_v_kadre(poziciya_kamery, lodka, bazis_kamery):
            continue
        
        zakryta = False
        for obrabotannaya_gora in obrabotannye_gori:
            if obrabotannaya_gora is None:
                continue
            if gora_zakryvaet_lodku(poziciya_kamery, lodka, obrabotannaya_gora):
                zakryta = True
                break
        
        if not zakryta:
            vidimye_lodki.append(nomer_lodki)
    
    return vidimye_lodki

def poluchit_poziciyu_kamery_na_progresse(nachalnaya_tochka, konechnaya_tochka, vektor_puti, progress):
    nachalo_x, nachalo_y, nachalo_z = nachalnaya_tochka
    put_dx, put_dy, put_dz = vektor_puti
    return (nachalo_x + put_dx * progress,
            nachalo_y + put_dy * progress,
            nachalo_z + put_dz * progress)

def poluchit_vidimye_lodki_na_progresse(nachalnaya_tochka, konechnaya_tochka, vektor_puti, progress,
                                   lodki, obrabotannye_gori, bazis_kamery):
    poziciya_kamery = poluchit_poziciyu_kamery_na_progresse(nachalnaya_tochka, konechnaya_tochka, vektor_puti, progress)
    return reshit_staticheskiy_sluchay(poziciya_kamery, lodki, obrabotannye_gori, bazis_kamery)

schetchik_vychisleniy = 0
kesh = {}
luchshee_kolichestvo_vidimyh = -1
luchshiy_progress = 0.0
luchshie_vidimye_lodki = []

def poluchit_vidimye_s_keshem(nachalnaya_tochka, konechnaya_tochka, vektor_puti, progress,
                               lodki, obrabotannye_gori, bazis_kamery):
    global schetchik_vychisleniy, kesh
    
    if progress in kesh:
        return kesh[progress]
    
    if schetchik_vychisleniy >= 1500:
        return 0, []
    
    vidimye_lodki = poluchit_vidimye_lodki_na_progresse(
        nachalnaya_tochka, konechnaya_tochka, vektor_puti, progress,
        lodki, obrabotannye_gori, bazis_kamery)
    
    kolichestvo = len(vidimye_lodki)
    kesh[progress] = (kolichestvo, vidimye_lodki)
    schetchik_vychisleniy += 1
    return kolichestvo, vidimye_lodki

def obnovit_luchshiy_rezultat(progress, kolichestvo_vidimyh, vidimye_lodki):
    global luchshee_kolichestvo_vidimyh, luchshiy_progress, luchshie_vidimye_lodki
    
    if kolichestvo_vidimyh > luchshee_kolichestvo_vidimyh or \
       (kolichestvo_vidimyh == luchshee_kolichestvo_vidimyh and progress < luchshiy_progress):
        luchshee_kolichestvo_vidimyh = kolichestvo_vidimyh
        luchshiy_progress = progress
        luchshie_vidimye_lodki = vidimye_lodki

def poisk_optimalnogo_momenta(nachalnaya_tochka, konechnaya_tochka, vektor_puti,
                              lodki, obrabotannye_gori, bazis_kamery,
                              levyy_progress, levoe_kolichestvo, levye_lodki,
                              pravyy_progress, pravoe_kolichestvo, pravye_lodki, glubina):
    global schetchik_vychisleniy
    
    if glubina >= 16 or schetchik_vychisleniy >= 1500:
        obnovit_luchshiy_rezultat(levyy_progress, levoe_kolichestvo, levye_lodki)
        obnovit_luchshiy_rezultat(pravyy_progress, pravoe_kolichestvo, pravye_lodki)
        return
    
    lokalnyy_maksimum = max(levoe_kolichestvo, pravoe_kolichestvo)
    if lokalnyy_maksimum < luchshee_kolichestvo_vidimyh:
        obnovit_luchshiy_rezultat(levyy_progress, levoe_kolichestvo, levye_lodki)
        obnovit_luchshiy_rezultat(pravyy_progress, pravoe_kolichestvo, pravye_lodki)
        return
    
    sredniy_progress = 0.5 * (levyy_progress + pravyy_progress)
    srednee_kolichestvo, srednie_lodki = poluchit_vidimye_s_keshem(
        nachalnaya_tochka, konechnaya_tochka, vektor_puti, sredniy_progress,
        lodki, obrabotannye_gori, bazis_kamery)
    obnovit_luchshiy_rezultat(sredniy_progress, srednee_kolichestvo, srednie_lodki)
    
    if max(levoe_kolichestvo, srednee_kolichestvo, pravoe_kolichestvo) < luchshee_kolichestvo_vidimyh:
        return
    
    poisk_optimalnogo_momenta(nachalnaya_tochka, konechnaya_tochka, vektor_puti,
                              lodki, obrabotannye_gori, bazis_kamery,
                              levyy_progress, levoe_kolichestvo, levye_lodki,
                              sredniy_progress, srednee_kolichestvo, srednie_lodki, glubina + 1)
    
    if schetchik_vychisleniy >= 1500:
        return
    
    poisk_optimalnogo_momenta(nachalnaya_tochka, konechnaya_tochka, vektor_puti,
                              lodki, obrabotannye_gori, bazis_kamery,
                              sredniy_progress, srednee_kolichestvo, srednie_lodki,
                              pravyy_progress, pravoe_kolichestvo, pravye_lodki, glubina + 1)

def nayti_optimalnyy(nachalnaya_tochka, konechnaya_tochka, vektor_puti, polnoe_vremya,
                     lodki, obrabotannye_gori, bazis_kamery):
    global schetchik_vychisleniy, kesh, luchshee_kolichestvo_vidimyh, luchshiy_progress, luchshie_vidimye_lodki
    
    schetchik_vychisleniy = 0
    kesh = {}
    luchshee_kolichestvo_vidimyh = -1
    luchshiy_progress = 0.0
    luchshie_vidimye_lodki = []
    
    levoe_kolichestvo, levye_lodki = poluchit_vidimye_s_keshem(
        nachalnaya_tochka, konechnaya_tochka, vektor_puti, 0.0,
        lodki, obrabotannye_gori, bazis_kamery)
    pravoe_kolichestvo, pravye_lodki = poluchit_vidimye_s_keshem(
        nachalnaya_tochka, konechnaya_tochka, vektor_puti, 1.0,
        lodki, obrabotannye_gori, bazis_kamery)
    
    obnovit_luchshiy_rezultat(0.0, levoe_kolichestvo, levye_lodki)
    obnovit_luchshiy_rezultat(1.0, pravoe_kolichestvo, pravye_lodki)
    
    poisk_optimalnogo_momenta(nachalnaya_tochka, konechnaya_tochka, vektor_puti,
                              lodki, obrabotannye_gori, bazis_kamery,
                              0.0, levoe_kolichestvo, levye_lodki,
                              1.0, pravoe_kolichestvo, pravye_lodki, 0)
    
    return luchshiy_progress * polnoe_vremya, luchshie_vidimye_lodki

def otformatirovat_vremya(vremya):
    if abs(vremya - round(vremya)) < 1e-9:
        return str(int(round(vremya)))
    else:
        return f"{vremya:.5f}"

def main():
    nachalnaya_tochka, konechnaya_tochka, skorost, napravlenie_kamery, lodki, gori = vvod_dannyh()
    
    bazis_kamery = postroit_bazis_kamery(napravlenie_kamery)
    if bazis_kamery is None:
        print("0")
        print(0)
        return
    
    obrabotannye_gori = predobrabotka_gor(gori)
    
    nachalo_x, nachalo_y, nachalo_z = nachalnaya_tochka
    konec_x, konec_y, konec_z = konechnaya_tochka
    
    put_dx = konec_x - nachalo_x
    put_dy = konec_y - nachalo_y
    put_dz = konec_z - nachalo_z
    dlina_puti = math.sqrt(put_dx * put_dx + put_dy * put_dy + put_dz * put_dz)
    
    if dlina_puti < 1e-9 or skorost <= 0:
        vidimye_lodki = reshit_staticheskiy_sluchay(nachalnaya_tochka, lodki, obrabotannye_gori, bazis_kamery)
        print("0")
        print(len(vidimye_lodki))
        for nomer_lodki in vidimye_lodki:
            print(nomer_lodki)
        return
    
    polnoe_vremya = dlina_puti / skorost
    vektor_puti = (put_dx, put_dy, put_dz)
    
    optimalnoe_vremya, vidimye_lodki = nayti_optimalnyy(nachalnaya_tochka, konechnaya_tochka, vektor_puti, polnoe_vremya,
                                 lodki, obrabotannye_gori, bazis_kamery)
    
    print(otformatirovat_vremya(optimalnoe_vremya))
    print(len(vidimye_lodki))
    for nomer_lodki in vidimye_lodki:
        print(nomer_lodki)

if __name__ == "__main__":
    main()
