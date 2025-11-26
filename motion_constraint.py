class MotionConstraint:
    def __init__(self, v0, v, a, j, s, t = 0):
        self.v0 = v0
        self.v  = v
        self.a  = a
        self.j  = j
        self.s  = s
        self.t  = t

    def update_a(self, a):
        self.a_prev = self.a
        self.a = a

    def restore_a(self):
        self.a = self.a_prev

    def update_v(self, v):
        self.v_prev = self.v
        self.v = v        

    def restore_v(self):
        self.v = self.v_prev
