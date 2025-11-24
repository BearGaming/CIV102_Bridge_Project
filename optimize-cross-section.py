import numpy as np
from math import sin, cos, tan, pi, radians, degrees
import matplotlib.pyplot as plt
from matplotlib.patches import Arc

def find_centroid(shapes):
    x = np.array([shape.centroid[0] for shape in shapes])
    y = np.array([shape.centroid[1] for shape in shapes])
    A = np.array([shape.area for shape in shapes])

    centroid = np.array([np.dot(A, x), np.dot(A, y)]) / np.sum(A)

    return centroid

class Quad:
    '''
    Pass position of the bottom left corner as a tuple
    Pass all lengths in mm
    Pass angles as degrees

    Defines a parallelogram of base (b), heigh (h) and interior angle of the bottom left corner (angle)
    '''
    def __init__(self, position, b, h, angle):
        # Defined as bottom left corner
        self.position = np.array(position)
        self.b = b
        self.h = h
        # Angle between slant and base
        self.angle = radians(angle)

        # Area and second moment of area
        self.area = b*h
        self.I = (b*(h**3))/12

        # Individual centroid
        self.slant = np.array([(self.h / tan(self.angle)), self.h])
        self.centroid = self.position + 0.5*(self.slant + np.array([self.b,0]))
        self.perimeter = 2*np.linalg.norm(self.slant) + 2*self.b

        # Coords
        # A list of the 4 vertices of the quadrilateral starting from the bottom left and going clockwise. 
        # The 5th entry is the first point to make drawing the lines between points easier
        self.coords = np.array([
                    self.position,
                    self.position + self.slant,
                    self.position + self.slant + np.array([self.b,0]),
                    self.position + np.array([self.b,0]),
                    self.position,
                ])

    def __str__(self):
        return f"y_b: {self.coords[0][1]}, y_t: {self.coords[1][1]}"

class CrossSection:
    def __init__(self, shapes):
        self.shapes = shapes
        self.centroid = find_centroid(shapes)
        self.d = []

        self.I = 0
        self.A = 0
        for shape in shapes:
            self.I += shape.I
            self.d.append(shape.centroid - self.centroid)
            self.I += (shape.area) * (((self.d[-1])[1])**2)
            
            self.A += (shape.area)
        self.perimeter = sum([shape.perimeter for shape in shapes])
        self.top = max([shape.coords[1][1] for shape in shapes])
    
    def __str__(self):
        desc = ""
        desc += "Height:" + str(self.top)
        desc += "\nCentroid:" + str(self.centroid) + "\nI: " + str(self.I)
        return desc

    def Q(self, depth_of_interest):
        Q = 0
        for shape in self.shapes:
            top = shape.coords[1][1]
            bottom = shape.coords[0][1]
            centroid = shape.centroid[1]

            area = 0
            d = 0

            if depth_of_interest > top:
                area = shape.area
                d = abs(self.centroid[1] - centroid)
            else:
                if depth_of_interest > bottom:
                    height = depth_of_interest - bottom
                    centroid = bottom + (height / 2)
                    area = height * shape.b
                    d = abs(self.centroid[1] - centroid)
            print(area, d, shape.b, "|", top, bottom, depth_of_interest)
            Q += area * d

        return Q
                
    def flexural_stress(self, M):
        y_top = self.top - self.centroid[1]
        y_bottom = self.centroid[1]

        fc = (M*y_top) / self.I
        ft = (M*y_bottom) / self.I

        return fc, ft


    def show(self, dist_to_centroid=False, title="", slant=False):
        fig = plt.figure()
        ax = fig.add_subplot()
        for shape in self.shapes:
            x = [coord[0] for coord in shape.coords]
            y = [coord[1] for coord in shape.coords]
            plt.plot(x, y, color="b")

            plt.plot(shape.centroid[0], shape.centroid[1], "rx")
        
        if dist_to_centroid:
            for d in self.d:
                value = d[1]
                d = d + np.array([0,self.centroid[1]])
                x = [d[0], d[0]]
                y = [self.centroid[1], d[1]]
                plt.plot(x, y, color="g")

        if slant:
            pos = slant.coords[2]
            angle = degrees(slant.angle)
            a = Arc(pos, 10, 10, angle=0, theta1=-angle, theta2=0, color='black', linewidth=1)
            angle = 180 - angle
            plt.text(pos[0] + 10, pos[1] - 10, str(round(angle, 3)) + '°')
            ax.add_patch(a)

        plt.plot(self.centroid[0], self.centroid[1], "rx", label=f"ȳ = {self.centroid[1]}")
        plt.hlines(self.centroid[1],xmin=-60,xmax=60, color='red', linestyle='--', linewidth=1)
        plt.legend(loc='lower right')
        plt.ylim(0, 120)
        plt.title(title)
        plt.axis('equal')
        plt.grid()
        plt.savefig(title)
        plt.show()

def build_trapezoid(angle, max_height, thickness, tab_width = 30, base=60, decks = 1):
    #assert base >= 60
    
    bottom = Quad((-base/2,0), base, thickness, 90)
    # Slanted sides
    left_slant = Quad((-base/2 - thickness, thickness), thickness, (max_height - ((2+decks)*thickness)), 90 + angle)
    right_slant = Quad((base/2, thickness), thickness, (max_height - ((2+decks)*thickness)), 90 - angle)

    # Glue tabs
    left_tab = Quad((left_slant.coords[1][0], left_slant.coords[1][1]), tab_width, thickness, 90)
    right_tab = Quad((right_slant.coords[1][0] - (tab_width - thickness), right_slant.coords[1][1]), tab_width, thickness, 90)

    U = [bottom, left_slant, right_slant, left_tab, right_tab]
    # Top
    top_width = (right_slant.coords[2] - left_slant.coords[1])[0]
    tops = [Quad((left_tab.coords[1][0], left_tab.coords[1][1]), top_width, thickness, 90)]
    for n in range(decks-1):
        tops.append(Quad((tops[n].coords[1][0],tops[n].coords[1][1]), top_width, thickness, 90))
    U.extend(tops)
    X = CrossSection(U)

    perimeter = bottom.b + top_width + 2*np.linalg.norm(left_slant.slant) + 2*tab_width
    return X, top_width, perimeter

def build_design_0():
    base = 80
    angle = 90
    max_height = 75+1.27
    thickness = 1.27
    decks = 1
    tab_width = 6.27
    
    bottom = Quad((-base/2,0), base, thickness, 90)
    # Slanted sides
    left_slant = Quad((-base/2 - thickness, thickness), thickness, (max_height - ((2+decks)*thickness)), 90+(90-angle))
    right_slant = Quad((base/2, thickness), thickness, (max_height - ((2+decks)*thickness)), angle)
    print(np.linalg.norm(left_slant.slant))

    # Glue tabs
    left_tab = Quad((left_slant.coords[1][0], left_slant.coords[1][1]), tab_width, thickness, 90)
    right_tab = Quad((right_slant.coords[1][0] - (tab_width - thickness), right_slant.coords[1][1]), tab_width, thickness, 90)

    U = [bottom, left_slant, right_slant, left_tab, right_tab]
    # Top
    top_width = 100
    tops = [Quad((left_tab.coords[1][0]-10, left_tab.coords[1][1]), top_width, thickness, 90)]
    U.extend(tops)
    X = CrossSection(U)

    perimeter = bottom.b + top_width + 2*np.linalg.norm(left_slant.slant) + 2*tab_width
    return X, perimeter, left_slant

def build_design_1():
    base = 100
    angle = 90
    thickness = 1.27
    max_height = 120
    decks = 2
    tab_width = 6.27
    
    bottom = Quad((-base/2,0), base, thickness, 90)
    # Slanted sides
    left_slant = Quad((-base/2, thickness), thickness, (max_height - ((2+decks)*thickness)), 90+(90-angle))
    right_slant = Quad((base/2 - thickness, thickness), thickness, (max_height - ((2+decks)*thickness)), angle)

    # Glue tabs
    left_tab = Quad((left_slant.coords[1][0], left_slant.coords[1][1]), tab_width, thickness, 90)
    right_tab = Quad((right_slant.coords[1][0] - (tab_width - thickness), right_slant.coords[1][1]), tab_width, thickness, 90)
    top_width = 100
    top1 = Quad((left_tab.coords[1][0], left_tab.coords[1][1]), top_width, thickness, 90)
    top2 = Quad((top1.coords[1][0], top1.coords[1][1]), top_width, thickness, 90)


    U = [bottom, left_slant, right_slant, left_tab, right_tab, top1, top2]
    # Top
    X = CrossSection(U)

    perimeter = bottom.b + top_width + 2*np.linalg.norm(left_slant.slant) + 2*tab_width
    return X, perimeter, left_slant

def build_design_2():
    base = 70
    angle = 82.644
    max_height = 120
    thickness = 1.27
    decks = 1
    tab_width = 30 + 1.27
    
    bottom = Quad((-base/2,0), base, thickness, 90)
    # Slanted sides
    left_slant = Quad((-base/2 - thickness, thickness), thickness, (max_height - ((2+decks)*thickness)), 90+(90-angle))
    right_slant = Quad((base/2, thickness), thickness, (max_height - ((2+decks)*thickness)), angle)
    print(np.linalg.norm(left_slant.slant))

    # Glue tabs
    left_tab = Quad((left_slant.coords[1][0], left_slant.coords[1][1]), tab_width, thickness, 90)
    right_tab = Quad((right_slant.coords[1][0] - (tab_width - thickness), right_slant.coords[1][1]), tab_width, thickness, 90)

    U = [bottom, left_slant, right_slant, left_tab, right_tab]
    # Top
    top_width = 100 + 2*1.27
    tops = [Quad((left_tab.coords[1][0], left_tab.coords[1][1]), top_width, thickness, 90)]
    for n in range(decks-1):
        tops.append(Quad((tops[n].coords[1][0],tops[n].coords[1][1]), top_width, thickness, 90))
    U.extend(tops)
    X = CrossSection(U)

    perimeter = bottom.b + top_width + 2*np.linalg.norm(left_slant.slant) + 2*tab_width
    return X, perimeter, left_slant

def build_design_3():
    base = 70
    angle = 81.136
    max_height = 100
    thickness = 1.27
    decks = 1
    tab_width = 30 + 1.27
    
    bottom = Quad((-base/2,0), base, thickness, 90)
    # Slanted sides
    left_slant = Quad((-base/2 - thickness, thickness), thickness, (max_height - ((2+decks)*thickness)), 90+(90-angle))
    right_slant = Quad((base/2, thickness), thickness, (max_height - ((2+decks)*thickness)), angle)
    print(np.linalg.norm(left_slant.slant))

    # Glue tabs
    left_tab = Quad((left_slant.coords[1][0], left_slant.coords[1][1]), tab_width, thickness, 90)
    right_tab = Quad((right_slant.coords[1][0] - (tab_width - thickness), right_slant.coords[1][1]), tab_width, thickness, 90)

    U = [bottom, left_slant, right_slant, left_tab, right_tab]
    # Top
    top_width = 100 + 2*1.27
    tops = [Quad((left_tab.coords[1][0], left_tab.coords[1][1]), top_width, thickness, 90)]
    for n in range(decks-1):
        tops.append(Quad((tops[n].coords[1][0],tops[n].coords[1][1]), top_width, thickness, 90))
    U.extend(tops)
    X = CrossSection(U)
    return X, 0, left_slant

def build_design_4():
    base = 70
    angle = 81.136
    max_height = 100
    thickness = 1.27
    decks = 1
    tab_width = 30 + 1.27
    
    bottom = Quad((-base/2,0), base, thickness, 90)
    # Slanted sides
    left_slant = Quad((-base/2 - thickness, thickness), thickness, (max_height - ((2+decks)*thickness)), 90+(90-angle))
    right_slant = Quad((base/2, thickness), thickness, (max_height - ((2+decks)*thickness)), angle)
    print(np.linalg.norm(left_slant.slant))

    # Glue tabs
    left_tab = Quad((left_slant.coords[1][0], left_slant.coords[1][1]), tab_width, thickness, 90)
    right_tab = Quad((right_slant.coords[1][0] - (tab_width - thickness), right_slant.coords[1][1]), tab_width, thickness, 90)

    U = [bottom, left_slant, right_slant, left_tab, right_tab]
    # Top
    top_width = 100 + 2*1.27
    tops = [Quad((left_tab.coords[1][0], left_tab.coords[1][1]), top_width, thickness, 90)]
    tops.append(Quad((left_tab.coords[3][0], left_tab.coords[3][1]), 40, thickness, 90))
    for n in range(decks-1):
        tops.append(Quad((tops[n].coords[1][0],tops[n].coords[1][1]), top_width, thickness, 90))
    U.extend(tops)
    X = CrossSection(U)


    perimeter = bottom.b + top_width + 2*np.linalg.norm(left_slant.slant) + 2*tab_width
    return X, perimeter, left_slant

def build_design_5a():
    base = 70
    angle = 81.13666
    max_height = 100
    thickness = 1.27
    tab_width = 50 + thickness
    
    bottom = Quad((-base/2,0), base, thickness, 90)
    # Slanted sides
    left_slant = Quad((-base/2 - thickness, thickness), thickness, (max_height - 3*thickness), 90+(90-angle))
    right_slant = Quad((base/2, thickness), thickness, (max_height - 3*thickness), angle)
    print(np.linalg.norm(left_slant.slant))
    
    left_tab = Quad((left_slant.coords[1][0], left_slant.coords[1][1]), tab_width, thickness, 90)
    right_tab = Quad((left_tab.coords[3][0], right_slant.coords[1][1]), tab_width, thickness, 90)

    bottom_width = 60
    r = Quad((left_slant.coords[1][0] + 21.27 , left_slant.coords[1][1] - thickness), bottom_width, thickness, 90)
    top_width = 100 + 2*thickness
    top2 = Quad((left_tab.coords[1][0], left_tab.coords[1][1]), top_width, thickness, 90)

    U = [bottom, left_slant, right_slant, r, left_tab, right_tab, top2]
    X = CrossSection(U)


    perimeter = bottom.b + top_width + 2*np.linalg.norm(left_slant.slant)
    return X, perimeter, left_slant

def build_design_5b():
    base = 70
    angle = 81.13666
    max_height = 100
    thickness = 1.27
    tab_width = 50 + thickness
    
    bottom = Quad((-base/2,0), base, thickness, 90)
    # Slanted sides
    left_slant = Quad((-base/2 - thickness, thickness), thickness, (max_height - 3*thickness), 90+(90-angle))
    right_slant = Quad((base/2, thickness), thickness, (max_height - 3*thickness), angle)
    print(np.linalg.norm(left_slant.slant))
    
    left_tab = Quad((left_slant.coords[1][0], left_slant.coords[1][1]), tab_width, thickness, 90)
    right_tab = Quad((left_tab.coords[3][0], right_slant.coords[1][1]), tab_width, thickness, 90)

    bottom_width = 90
    r = Quad((left_slant.coords[1][0] + 6.27 , left_slant.coords[1][1] - thickness), bottom_width, thickness, 90)
    top_width = 100 + 2*thickness
    top2 = Quad((left_tab.coords[1][0], left_tab.coords[1][1]), top_width, thickness, 90)

    U = [bottom, left_slant, right_slant, r, left_tab, right_tab, top2]
    X = CrossSection(U)


    perimeter = bottom.b + top_width + 2*np.linalg.norm(left_slant.slant)
    return X, perimeter, left_slant

def optimize_trapezoid():
    # Iterate through heights in multiples of 20
    for h in range(20, 201, 20):
        print("#"*20)
        print(f"Height = {h}mm")
        # Iterate through angles starting from shallow to steep
        for i in range(89,1,-1):
            X, width, perimeter = build_trapezoid(i, h, 1.27)
            if width < 100:
                print("deck width too small")
                break

            load_ratio = X.centroid[1] / (X.top - X.centroid[1])
            if load_ratio < 5:
                print("Angle:", 90 - i, "Load Ratio:", load_ratio, "Deck width:", width, "I:", X.I)
        X.show(title="Height = " + str(h))

def show_designs():
    X, p, a = build_design_1()
    print(X)
    X.show(title="Design 1", slant = a)
    
    X, p, a = build_design_2()
    print(X)
    X.show(title="Design 2", slant = a)
    X, p, a = build_design_3()
    print(X)
    X.show(title="Design 3", slant = a)
    X, p, a = build_design_4()
    print(X)
    X.show(title="Design 4", slant = a)
    X, p, a = build_design_5b()
    print(X)
    X.show(title="Design 5", slant = a)

if __name__ == "__main__":
    optimize_trapezoid()