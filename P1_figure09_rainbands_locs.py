from wrf import CoordPair
def get_rainbands_locs():
 start_point1 = CoordPair(lat=29.7, lon=-95.6) # for 2512
 end_point1 = CoordPair(lat=30.22, lon=-94.5)  # for 2512
 rainband1 = {"start": start_point1, "end": end_point1}
 
 start_point2 = CoordPair(lat=28.74, lon=-96.09) # for 2512
 end_point2 = CoordPair(lat=29.37, lon=-94.98)  # for 2512
 rainband2 = {"start": start_point2, "end": end_point2}
 return {'rainband1':rainband1, 'rainband2':rainband2} 




