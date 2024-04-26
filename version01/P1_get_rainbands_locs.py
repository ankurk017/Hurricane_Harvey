from wrf import CoordPair
def get_rainbands_locs():
 start_point1 = CoordPair(lat=29.7, lon=-95.6) # for 2512
 end_point1 = CoordPair(lat=30.22, lon=-94.5)  # for 2512
 rainband1 = {"start": start_point1, "end": end_point1}
 
 start_point2 = CoordPair(lat=28.74, lon=-96.09) # for 2512
 end_point2 = CoordPair(lat=29.37, lon=-94.98)  # for 2512
 rainband2 = {"start": start_point2, "end": end_point2}
 return {'rainband1':rainband1, 'rainband2':rainband2} 


def get_rainbands_locs_updated():
 start_point1 = CoordPair(lat=29.87, lon=-95.5) # for 2512
 end_point1 = CoordPair(lat=30.03, lon=-94.31)  # for 2512

 start_point1 = CoordPair(lat=29.8, lon=-95.89) # for 2512
 end_point1 = CoordPair(lat=30.15, lon=-93.92)  # for 2512
 rainband1 = {"start": start_point1, "end": end_point1}

 start_point2 = CoordPair(lat=29.02, lon=-95.94) # for 2512
 end_point2 = CoordPair(lat=29.02, lon=-94.55)  # for 2512
 rainband4 = {"start": start_point2, "end": end_point2}

 start_point1 = CoordPair(lat=29.9, lon=-95.87) # for 2512
 end_point1 = CoordPair(lat=29.9, lon=-94.22)  # for 2512
 rainband3 = {"start": start_point1, "end": end_point1}

 start_point2 = CoordPair(lat=28.64, lon=-96.01) # for 2512
 end_point2 = CoordPair(lat=29.44, lon=-94.71)  # for 2512
 rainband2 = {"start": start_point2, "end": end_point2}


 return {'rainband1':rainband1, 'rainband2':rainband2, 'rainband3':rainband3, 'rainband4':rainband4}






