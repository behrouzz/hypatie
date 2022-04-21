class RAhms:
    def __init__(self, h=0, m=0, s=0):
        self.h = h
        self.m = m
        self.s = s
        self.hours = h + m/60 + s/3600
        self.minutes = h*60 + m + s/60
        self.seconds = h*3600 + m*60 + s
        self.deg = 15 * self.hours
        self.rad = self.deg * (3.141592653589793/180)
        

class DECdms:
    def __init__(self, sign='+', d=0, m=0, s=0):
        self.sign = sign
        sgn = -1 if sign=='-' else 1
        self.d = d
        self.m = m
        self.s = s
        self.degrees = (d + m/60 + s/3600) * sgn
        self.minutes = (d*60 + m + s/60) * sgn
        self.seconds = (d*3600 + m*60 + s) * sgn
        self.deg = self.degrees
        self.rad = self.deg * (3.141592653589793/180)


class RAdeg:
    def __init__(self, ra):
        self.rad = ra * (3.141592653589793/180)
        self.h = int(ra/15)
        self.m = int(((ra/15)-self.h)*60)
        self.s = ((((ra/15)-self.h)*60)-self.m)*60
        self.str = str(self.h)+'h' + str(self.m)+'m' + str(self.s)+'s'
        

class DECdeg:
    def __init__(self, dec):
        self.rad = dec * (3.141592653589793/180)
        self.sign = '-' if dec<0 else '+'
        dec = abs(dec)
        self.d = int(dec)
        self.m = abs(int((dec-self.d)*60))
        self.s = (abs((dec-self.d)*60)-self.m)*60
        self.str = self.sign + str(self.d)+'d' + \
                   str(self.m)+'m' + str(self.s)+'s'
