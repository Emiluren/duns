use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

#[derive(Clone, Copy, Debug)]
pub struct Vec3 {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

#[derive(Clone, Debug)]
pub struct Mat3 {
    pub data: [f32; 9]
}

#[derive(Clone, Debug)]
pub struct Mat4 {
    pub data: [f32; 12]
}

#[derive(Clone, Copy, Debug)]
pub struct Quaternion {
    pub r: f32,
    pub i: f32,
    pub j: f32,
    pub k: f32,
}

pub fn vec3(x: f32, y: f32, z: f32) -> Vec3 {
    Vec3 { x, y, z }
}

impl Add for Vec3 {
    type Output = Vec3;

    fn add(self, other: Vec3) -> Self {
        vec3(
            self.x + other.x,
            self.y + other.y,
            self.z + other.z
        )
    }
}

impl AddAssign for Vec3 {
    fn add_assign(&mut self, other: Vec3) {
        *self = *self + other;
    }
}

impl Sub for Vec3 {
    type Output = Vec3;

    fn sub(self, other: Vec3) -> Self {
        vec3(
            self.x - other.x,
            self.y - other.y,
            self.z - other.z,
        )
    }
}

impl SubAssign for Vec3 {
    fn sub_assign(&mut self, other: Vec3) {
        *self = *self - other;
    }
}

impl Neg for Vec3 {
    type Output = Vec3;

    fn neg(self) -> Self {
        vec3(
            -self.x,
            -self.y,
            -self.z,
        )
    }
}

impl Mul<f32> for Vec3 {
    type Output = Vec3;

    fn mul(self, scalar: f32) -> Self {
        vec3(
            self.x * scalar,
            self.y * scalar,
            self.z * scalar,
        )
    }
}

impl MulAssign<f32> for Vec3 {
    fn mul_assign(&mut self, scalar: f32) {
        *self = *self * scalar;
    }
}

impl Mul<Vec3> for Vec3 {
    type Output = Vec3;

    fn mul(self, other: Self) -> Self {
        vec3(
            self.x * other.x,
            self.y * other.y,
            self.z * other.z,
        )
    }
}

impl MulAssign<Vec3> for Vec3 {
    fn mul_assign(&mut self, other: Vec3) {
        *self = *self * other;
    }
}

impl Div<f32> for Vec3 {
    type Output = Vec3;

    fn div(self, scalar: f32) -> Self {
        vec3(
            self.x / scalar,
            self.y / scalar,
            self.z / scalar,
        )
    }
}

impl DivAssign<f32> for Vec3 {
    fn div_assign(&mut self, scalar: f32) {
        *self = *self / scalar
    }
}

impl Vec3 {
    pub fn norm(self) -> f32 {
        (self.x.powi(2) + self.y.powi(2) + self.z.powi(2)).sqrt()
    }

    pub fn normalize(self) -> Vec3 {
        self / self.norm()
    }

    pub fn distance_to(self, other: Self) -> f32 {
        (self - other).norm()
    }

    pub fn dot(self, other: Self) -> f32 {
        self.x * other.x + self.y * other.y
    }

    pub fn cross(self, other: Self) -> Vec3 {
        vec3(
            self.y*other.z - self.z*other.y,
            self.z*other.x - self.x*other.z,
            self.x*other.y - self.y*other.x,
        )
    }

    pub fn to_world(self, transform: Mat4) -> Vec3 {
        transform * self
    }

    pub fn to_local(self, transform: Mat4) -> Vec3 {
        transform.transform_inverse(self)
    }

    pub fn to_world_dir(self, transform: Mat4) -> Vec3 {
        transform.transform_direction(self)
    }

    pub fn to_local_dir(self, transform: Mat4) -> Vec3 {
        transform.transform_inverse_direction(self)
    }
}

impl Mat3 {
    pub fn from_orientation(q: Quaternion) -> Mat3 {
        Mat3 { data: [
            1.0 - (2.0*q.j*q.j + 2.0*q.k*q.k),
            2.0*q.i*q.j + 2.0*q.k*q.r,
            2.0*q.i*q.k - 2.0*q.j*q.r,
            2.0*q.i*q.j - 2.0*q.k*q.r,
            1.0 - (2.0*q.i*q.i + 2.0*q.k*q.k),
            2.0*q.j*q.k + 2.0*q.i*q.r,
            2.0*q.i*q.k + 2.0*q.j*q.r,
            2.0*q.j*q.k - 2.0*q.i*q.r,
            1.0 - (2.0*q.i*q.i + 2.0*q.j*q.j),
        ]}
    }

    pub fn row0(&self) -> Vec3 {
        vec3(
            self.data[0],
            self.data[1],
            self.data[2],
        )
    }

    pub fn row1(&self) -> Vec3 {
        vec3(
            self.data[3],
            self.data[4],
            self.data[5],
        )
    }

    pub fn row2(&self) -> Vec3 {
        vec3(
            self.data[6],
            self.data[7],
            self.data[8],
        )
    }

    pub fn col0(&self) -> Vec3 {
        vec3(
            self.data[0],
            self.data[3],
            self.data[6],
        )
    }

    pub fn col1(&self) -> Vec3 {
        vec3(
            self.data[1],
            self.data[4],
            self.data[7],
        )
    }

    pub fn col2(&self) -> Vec3 {
        vec3(
            self.data[2],
            self.data[5],
            self.data[8],
        )
    }

    pub fn invert(&mut self) {
        let t4 = self.data[0]*self.data[4];
        let t6 = self.data[0]*self.data[5];
        let t8 = self.data[1]*self.data[3];
        let t10 = self.data[2]*self.data[3];
        let t12 = self.data[1]*self.data[6];
        let t14 = self.data[2]*self.data[6];

        // Calculate the determinant.
        let t16 = t4*self.data[8] - t6*self.data[7] - t8*self.data[8] +
            t10*self.data[7] + t12*self.data[5] - t14*self.data[4];

        // Make sure the determinant is non-zero.
        if t16 == 0.0 {
            return;
        }
        let t17 = 1.0 / t16;

        self.data[0] = (self.data[4]*self.data[8]-self.data[5]*self.data[7])*t17;
        self.data[1] = -(self.data[1]*self.data[8]-self.data[2]*self.data[7])*t17;
        self.data[2] = (self.data[1]*self.data[5]-self.data[2]*self.data[4])*t17;
        self.data[3] = -(self.data[3]*self.data[8]-self.data[5]*self.data[6])*t17;
        self.data[4] = (self.data[0]*self.data[8]-t14)*t17;
        self.data[5] = -(t6-t10)*t17;
        self.data[6] = (self.data[3]*self.data[7]-self.data[4]*self.data[6])*t17;
        self.data[7] = -(self.data[0]*self.data[7]-t12)*t17;
        self.data[8] = (t4-t8)*t17;
    }

    pub fn inverse(&self) -> Mat3 {
        let mut res = self.clone();
        res.invert();
        res
    }

    pub fn set_transpose(&mut self) {
        self.data[0] = self.data[0];
        self.data[1] = self.data[3];
        self.data[2] = self.data[6];
        self.data[3] = self.data[1];
        self.data[4] = self.data[4];
        self.data[5] = self.data[7];
        self.data[6] = self.data[2];
        self.data[7] = self.data[5];
        self.data[8] = self.data[8];
    }

    pub fn transpose(&self) -> Mat3 {
        let mut res = self.clone();
        res.set_transpose();
        res
    }
}

impl Mul<Vec3> for Mat3 {
    type Output = Vec3;

    fn mul(self, rhs: Vec3) -> Vec3 {
        vec3(
            self.row0().dot(rhs),
            self.row1().dot(rhs),
            self.row2().dot(rhs),
        )
    }
}

impl Mul<Mat3> for Mat3 {
    type Output = Mat3;

    fn mul(self, rhs: Mat3) -> Mat3 {
        Mat3 { data: [
            self.row0().dot(rhs.col0()),
            self.row0().dot(rhs.col1()),
            self.row0().dot(rhs.col2()),

            self.row1().dot(rhs.col0()),
            self.row1().dot(rhs.col1()),
            self.row1().dot(rhs.col2()),

            self.row2().dot(rhs.col0()),
            self.row2().dot(rhs.col1()),
            self.row2().dot(rhs.col2()),
        ]}
    }
}

impl Mul<Vec3> for Mat4 {
    type Output = Vec3;

    fn mul(self, rhs: Vec3) -> Vec3 {
        vec3(
            rhs.x*self.data[0] + rhs.y*self.data[1] + rhs.z*self.data[2] + self.data[3],
            rhs.x*self.data[4] + rhs.y*self.data[5] + rhs.z*self.data[6] + self.data[7],
            rhs.x*self.data[8] + rhs.y*self.data[9] + rhs.z*self.data[10] + self.data[11],
        )
    }
}

impl Mat4 {
    pub fn from_orientation_and_pos(q: Quaternion, pos: Vec3) -> Mat4 {
        Mat4 { data: [
            1.0 - (2.0*q.j*q.j + 2.0*q.k*q.k),
            2.0*q.i*q.j + 2.0*q.k*q.r,
            2.0*q.i*q.k - 2.0*q.j*q.r,
            pos.x,
            2.0*q.i*q.j - 2.0*q.k*q.r,
            1.0 - (2.0*q.i*q.i + 2.0*q.k*q.k),
            2.0*q.j*q.k + 2.0*q.i*q.r,
            pos.y,
            2.0*q.i*q.k + 2.0*q.j*q.r,
            2.0*q.j*q.k - 2.0*q.i*q.r,
            1.0 - (2.0*q.i*q.i + 2.0*q.j*q.j),
            pos.z,
        ]}
    }

    pub fn get_determinant(&self) -> f32 {
        self.data[8]*self.data[5]*self.data[2]+
            self.data[4]*self.data[9]*self.data[2]+
            self.data[8]*self.data[1]*self.data[6]-
            self.data[0]*self.data[9]*self.data[6]-
            self.data[4]*self.data[1]*self.data[10]+
            self.data[0]*self.data[5]*self.data[10]
    }

    pub fn invert(&mut self) {
        let det = self.get_determinant();
        if det == 0.0 {
            return;
        }
        let inv_det = 1.0 / det;

        self.data[0] = (-self.data[9]*self.data[6]+self.data[5]*self.data[10])*inv_det;
        self.data[4] = (self.data[8]*self.data[6]-self.data[4]*self.data[10])*inv_det;
        self.data[8] = (-self.data[8]*self.data[5]+self.data[4]*self.data[9]* self.data[15])*inv_det;

        self.data[1] = (self.data[9]*self.data[2]-self.data[1]*self.data[10])*inv_det;
        self.data[5] = (-self.data[8]*self.data[2]+self.data[0]*self.data[10])*inv_det;
        self.data[9] = (self.data[8]*self.data[1]-self.data[0]*self.data[9]* self.data[15])*inv_det;

        self.data[2] = (-self.data[5]*self.data[2]+self.data[1]*self.data[6]* self.data[15])*inv_det;
        self.data[6] = (self.data[4]*self.data[2]-self.data[0]*self.data[6]* self.data[15])*inv_det;
        self.data[10] = (-self.data[4]*self.data[1]+self.data[0]*self.data[5]* self.data[15])*inv_det;

        self.data[3] = (
            self.data[9]*self.data[6]*self.data[3] -
                self.data[5]*self.data[10]*self.data[3] -
                self.data[9]*self.data[2]*self.data[7] +
                self.data[1]*self.data[10]*self.data[7] +
                self.data[5]*self.data[2]*self.data[11] -
                self.data[1]*self.data[6]*self.data[11]
        )*inv_det;
        self.data[7] = (
            -self.data[8]*self.data[6]*self.data[3] +
                self.data[4]*self.data[10]*self.data[3] +
                self.data[8]*self.data[2]*self.data[7] -
                self.data[0]*self.data[10]*self.data[7] -
                self.data[4]*self.data[2]*self.data[11] +
                self.data[0]*self.data[6]*self.data[11]
        )*inv_det;
        self.data[11] =(
            self.data[8]*self.data[5]*self.data[3]
                -self.data[4]*self.data[9]*self.data[3] -
                self.data[8]*self.data[1]*self.data[7] +
                self.data[0]*self.data[9]*self.data[7] +
                self.data[4]*self.data[1]*self.data[11] -
                self.data[0]*self.data[5]*self.data[11]
        )*inv_det;
    }

    pub fn inverse(&self) -> Mat4 {
        let mut res = self.clone();
        res.invert();
        res
    }

    pub fn transform_inverse(&self, v: Vec3) -> Vec3 {
        let tmp = v - vec3(self.data[3], self.data[7], self.data[11]);
        vec3(
            tmp.dot(vec3(self.data[0], self.data[4], self.data[8])),
            tmp.dot(vec3(self.data[1], self.data[5], self.data[9])),
            tmp.dot(vec3(self.data[2], self.data[6], self.data[10])),
        )
    }

    pub fn transform_direction(&self, v: Vec3) -> Vec3 {
        vec3(
            v.dot(vec3(self.data[0], self.data[1], self.data[2])),
            v.dot(vec3(self.data[4], self.data[5], self.data[6])),
            v.dot(vec3(self.data[8], self.data[9], self.data[10])),
        )
    }

    pub fn transform_inverse_direction(&self, v: Vec3) -> Vec3 {
        vec3(
            v.dot(vec3(self.data[0], self.data[4], self.data[8])),
            v.dot(vec3(self.data[1], self.data[5], self.data[9])),
            v.dot(vec3(self.data[2], self.data[6], self.data[10])),
        )
    }
}

impl Mul<Mat4> for Mat4 {
    type Output = Mat4;

    fn mul(self, rhs: Mat4) -> Mat4 {
        let l = &self.data;
        let r = &rhs.data;
        Mat4 { data: [
            l[0]*r[0] + l[1]*r[4] + l[2]*r[8],
            l[0]*r[1] + l[1]*r[5] + l[2]*r[9],
            l[0]*r[2] + l[1]*r[6] + l[2]*r[10],
            l[0]*r[3] + l[1]*r[7] + l[2]*r[11] + l[3],

            l[4]*r[0] + l[5]*r[4] + l[6]*r[8],
            l[4]*r[1] + l[5]*r[5] + l[6]*r[9],
            l[4]*r[2] + l[5]*r[6] + l[6]*r[10],
            l[4]*r[3] + l[5]*r[7] + l[6]*r[11] + l[7],

            l[8]*r[0] + l[9]*r[4] + l[10]*r[8],
            l[8]*r[1] + l[9]*r[5] + l[10]*r[9],
            l[8]*r[2] + l[9]*r[6] + l[10]*r[10],
            l[8]*r[3] + l[9]*r[7] + l[10]*r[11] + l[11],
        ]}
    }
}

impl Quaternion {
    pub fn normalize(&mut self) {
        let d2 = self.r*self.r + self.i*self.i + self.j*self.j + self.k*self.k;

        // Check for zero length quaternion and use the no-rotation quaternion for that case
        if d2 == 0.0 {
            self.r = 1.0;
            return;
        }

        let d = d2.sqrt();
        self.r /= d;
        self.i /= d;
        self.j /= d;
        self.k /= d;
    }

    pub fn rotate_by_vector(self, v: Vec3) -> Quaternion {
        self * Quaternion {
            r: 0.0,
            i: v.x,
            j: v.y,
            k: v.z,
        }
    }

    pub fn add_scaled_vector(&mut self, v: Vec3, scale: f32) {
        let q = Quaternion { r: 0.0, i: v.x, j: v.y, k: v.z } * scale;
        *self = q * (*self) * 0.5;
    }
}

impl Mul<Quaternion> for Quaternion {
    type Output = Quaternion;

    fn mul(self, rhs: Quaternion) -> Quaternion {
        Quaternion {
            r: self.r*rhs.r - self.i*rhs.i - self.j*rhs.j - self.k*rhs.k,
            i: self.r*rhs.i + self.i*rhs.r + self.j*rhs.k - self.k*rhs.j,
            j: self.r*rhs.j + self.j*rhs.r + self.k*rhs.i - self.i*rhs.k,
            k: self.r*rhs.k + self.k*rhs.r + self.i*rhs.j - self.j*rhs.i,
        }
    }
}

impl MulAssign<Quaternion> for Quaternion {
    fn mul_assign(&mut self, rhs: Quaternion) {
        *self = *self * rhs;
    }
}

impl Mul<f32> for Quaternion {
    type Output = Quaternion;

    fn mul(self, rhs: f32) -> Quaternion {
        Quaternion {
            r: self.r * rhs,
            i: self.i * rhs,
            j: self.j * rhs,
            k: self.k * rhs,
        }
    }
}
