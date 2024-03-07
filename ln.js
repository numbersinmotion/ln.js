let ln = {
	INF: Number.POSITIVE_INFINITY,
	EPS: Number.EPSILON,
	AXIS: { NONE: "AxisNone", X: "AxisX", Y: "AxisY", Z: "AxisZ" }
};

ln.Utils = {
	degToRad(degrees) { return degrees * Math.PI / 180; },
	radToDeg(radians) { return radians * 180 / Math.PI; },
	getMedian: function (items) {
		let n = items.length;
		if (n === 0) return 0;
		if (n % 2 === 1) return items[Math.ceil(n / 2)];
		return (items[n / 2] + items[n / 2 - 1]) / 2;
	}
}

ln.Scene = class {
	constructor() { [this.Shapes, this.Tree] = [[], null]; }
	compile() {
		for (let shape of this.Shapes) shape.Compile();
		if (this.Tree === null) this.Tree = ln.Tree.NewTree(this.Shapes);
	}
	addShape(shape) { this.Shapes.push(shape); }
	Intersect(ray) { return this.Tree.Intersect(ray); }
	Visible(eye, point) {
		let v = eye.Sub(point);
		let r = new ln.Ray(point, v.Normalize());
		let hit = this.Intersect(r);
		return hit.T >= v.Length();
	}
	Paths() {
		let paths = new ln.Paths();
		for (let shape of this.Shapes) {
			for (let path of shape.Paths().paths) paths.Add(path);
		}
		return paths;
	}
	Render(eye, center, up, width, height, fovy, near, far, step) {
		let [aspect, matrix] = [width / height, ln.Matrix.createLookAt(eye, center, up)];
		matrix = matrix.Perspective(fovy, aspect, near, far);
		return this.RenderWithMatrix(matrix, eye, width, height, step);
	}
	RenderWithMatrix(matrix, eye, width, height, step) {
		this.compile();
		let paths = this.Paths();
		if (step > 0) paths = paths.Chop(step);
		paths = paths.Filter(new ln.ClipFilter(matrix, eye, this));
		if (step > 0) paths = paths.Simplify(1e-6);
		matrix = ln.Matrix.createTranslate(new ln.Vector(1, 1, 0)).Scale(new ln.Vector(0.5 * width, 0.5 * height, 0));
		paths = paths.Transform(matrix);
		return paths;
	}
};

ln.Hit = class {
	constructor(shape, t) { [this.Shape, this.T] = [shape, t]; }
	static none() { return new ln.Hit(null, ln.INF); }
	// Ok() { return this.T < ln.INF; }
	min(h) { return this.T <= h.t ? this : h; }
	max(h) { return this.T > h.T ? this : h; }
};

ln.ClipFilter = class {
	constructor(matrix, eye, scene) {
		[this.Matrix, this.Eye, this.Scene] = [matrix, eye, scene];
		this.ClipBox = new ln.Box(new ln.Vector(-1, -1, -1), new ln.Vector(1, 1, 1));
	}
	Filter(v) {
		let w = this.Matrix.MulPositionW(v);
		if (!this.Scene.Visible(this.Eye, v)) return [w, false];
		if (!this.ClipBox.Contains(w)) return [w, false];
		return [w, true];
	}
};

ln.Path = class {
	constructor() { this.points = []; }
	static FromPoints(points) {
		let path = new ln.Path();
		path.points = points;
		return path;
	}
	Add(v) { this.points.push(v); }
	BoundingBox() {
		let box = new ln.Box(this.points[0], this.points[0]);
		for (let i = 1; i < this.points.length; i++) box = box.Extend(new ln.Box(this.points[i], this.points[i]));
		return box;
	}
	Transform(matrix) {
		let path = new ln.Path();
		for (let p of this.points) path.Add(matrix.MulPosition(p));
		return path;
	}
	Chop(step) {
		let path = new ln.Path();
		for (let i = 0; i < this.points.length - 1; i++) {
			let [a, b] = [this.points[i], this.points[i + 1]];
			let v = b.Sub(a);
			let l = v.Length();
			if (i === 0) path.Add(a);
			let d = step;
			for (let d = step; d < l; d += step) path.Add(a.Add(v.MulScalar(d / l)));
			path.Add(b);
		}
		return path;
	}
	Filter(f) {
		let [path, paths] = [new ln.Path(), new ln.Paths()];
		for (let p of this.points) {
			let [v, ok] = f.Filter(p);
			if (ok) path.Add(v);
			else {
				if (path.points.length > 1) paths.Add(path);
				path = new ln.Path();
			}
		}
		if (path.points.length > 1) paths.Add(path);
		return paths;
	}
	Simplify(threshold) {
		if (this.points.length < 3) return this;
		let [a, b, index, distance] = [this.points[0], this.points[this.points.length - 1], -1, 0];
		for (let i = 1; i < this.points.length - 1; i++) {
			let d = this.points[i].SegmentDistance(a, b);
			if (d > distance) [index, distance] = [i, d];
		}
		if (distance > threshold) {
			let r1 = ln.Path.FromPoints(this.points.slice(0, index + 1)).Simplify(threshold);
			let r2 = ln.Path.FromPoints(this.points.slice(index)).Simplify(threshold);
			return ln.Path.FromPoints([...r1.points.slice(0, r1.points.length - 1), ...r2.points]);
		} else {
			return ln.Path.FromPoints([a, b]);
		}
	}
}

ln.Paths = class {
	constructor() { this.paths = []; }
	Add(path) { this.paths.push(path); }
	BoundingBox() {
		let box = this.paths[0].BoundingBox();
		for (let i = 1; i < this.paths.legnth; i++) box = box.Extend(this.paths[i].BoundingBox());
		return box;
	}
	Transform(matrix) {
		let paths = new ln.Paths();
		for (let path of this.paths) paths.Add(path.Transform(matrix));
		return paths;
	}
	Chop(step) {
		let paths = new ln.Paths();
		for (let path of this.paths) paths.Add(path.Chop(step));
		return paths;
	}
	Filter(f) {
		let paths = new ln.Paths();
		for (let path of this.paths) {
			for (let filteredPath of path.Filter(f).paths) paths.Add(filteredPath);
		}
		return paths;
	}
	Simplify(threshold) {
		let paths = new ln.Paths();
		for (let path of this.paths) paths.Add(path.Simplify(threshold));
		return paths;
	}
}

ln.Tree = class {
	constructor(box, node) { [this.Box, this.Root] = [box, node]; }
	static NewTree(shapes) {
		let box = ln.Box.BoxForShapes(shapes);
		let node = ln.Node.NewNode(shapes);
		node.Split(0);
		return new ln.Tree(box, node);
	}
	Intersect(r) {
		let [tmin, tmax] = this.Box.Intersect(r);
		if (tmax < tmin || tmax <= 0) return ln.Hit.none();
		return this.Root.Intersect(r, tmin, tmax);
	}
};

ln.Node = class {
	constructor(axis, point, shapes, left, right) { [this.Axis, this.Point, this.Shapes, this.Left, this.Right] = [axis, point, shapes, left, right]; }
	static NewNode(shapes) { return new ln.Node(ln.AXIS.NONE, 0, shapes, null, null); }
	Intersect(r, tmin, tmax) {
		let tsplit, leftFirst;
		switch (this.Axis) {
			case ln.AXIS.NONE:
				return this.IntersectShapes(r);
				break;
			case ln.AXIS.X:
				tsplit = (this.Point - r.Origin.X) / r.Direction.X;
				leftFirst = (r.Origin.X < this.Point) || (r.Origin.X === this.Point && r.Direction.X <= 0);
				break;
			case ln.AXIS.Y:
				tsplit = (this.Point - r.Origin.Y) / r.Direction.Y;
				leftFirst = (r.Origin.Y < this.Point) || (r.Origin.Y === this.Point && r.Direction.Y <= 0);
				break;
			case ln.AXIS.Z:
				tsplit = (this.Point - r.Origin.Z) / r.Direction.Z;
				leftFirst = (r.Origin.Z < this.Point) || (r.Origin.Z === this.Point && r.Direction.Z <= 0);
				break;
		}
		let [first, second] = leftFirst ? [this.Left, this.Right] : [this.Right, this.Left];
		if (tsplit > tmax || tsplit <= 0) {
			return first.Intersect(r, tmin, tmax);
		} else if (tsplit < tmin) {
			return second.Intersect(r, tmin, tmax);
		} else {
			let h1 = first.Intersect(r, tmin, tsplit);
			if (h1.T <= tsplit) return h1;
			let h2 = second.Intersect(r, tsplit, Math.min(tmax, h1.T));
			if (h1.T <= h2.T) return h1;
			else return h2;
		}
	}
	IntersectShapes(r) {
		let hit = ln.Hit.none();
		for (let shape of this.Shapes) {
			let h = shape.Intersect(r);
			if (h.T < hit.T) hit = h;
		}
		return hit;
	}
	PartitionScore(axis, point) {
		let [left, right] = [0, 0];
		for (let shape of this.Shapes) {
			let box = shape.BoundingBox();
			let [l, r] = box.Partition(axis, point);
			if (l) left++;
			if (r) right++;
		}
		if (left >= right) return left;
		return right;
	}
	Partition(size, axis, point) {
		let [left, right] = [[], []];
		for (let shape of this.Shapes) {
			let box = shape.BoundingBox();
			let [l, r] = box.Partition(axis, point);
			if (l) left.push(shape);
			if (r) right.push(shape);
		}
		return [left, right];
	}
	Split(depth) {
		if (this.Shapes.length < 8) return;
		let [xs, ys, zs] = [[], [], []];
		for (let shape of this.Shapes) {
			let box = shape.BoundingBox();
			xs.push(box.Min.X, box.Max.X);
			ys.push(box.Min.Y, box.Max.Y);
			zs.push(box.Min.Z, box.Max.Z);
		}
		let comp = function (a, b) { return a - b; };
		xs.sort(comp);
		ys.sort(comp);
		zs.sort(comp);
		let [mx, my, mz] = [ln.Utils.getMedian(xs), ln.Utils.getMedian(ys), ln.Utils.getMedian(zs)];
		let [bestScore, bestAxis, bestPoint] = [Math.floor(this.Shapes.length * 0.85), ln.AXIS.NONE, 0];
		let sx = this.PartitionScore(ln.AXIS.X, mx);
		if (sx < bestScore) [bestScore, bestAxis, bestPoint] = [sx, ln.AXIS.X, mx];
		let sy = this.PartitionScore(ln.AXIS.Y, my);
		if (sy < bestScore) [bestScore, bestAxis, bestPoint] = [sy, ln.AXIS.Y, my];
		let sz = this.PartitionScore(ln.AXIS.Z, mz);
		if (sz < bestScore) [bestScore, bestAxis, bestPoint] = [sz, ln.AXIS.Z, mz];
		if (bestAxis === ln.AXIS.NONE) return;
		let [l, r] = this.Partition(bestScore, bestAxis, bestPoint);
		[this.Axis, this.Point, this.Left, this.Right] = [bestAxis, bestPoint, ln.Node.NewNode(l), ln.Node.NewNode(r)];
		this.Left.Split(depth + 1);
		this.Right.Split(depth + 1);
		this.Shapes = null;
	}
}

ln.Box = class {
	constructor(min, max) { [this.Min, this.Max] = [min, max]; }
	static BoxForShapes(shapes) {
		if (shapes.length === 0) return new ln.Box(new ln.Vector(0, 0, 0), new ln.Vector(0, 0, 0));
		let box = shapes[0].BoundingBox();
		for (let i = 1; i < shapes.length; i++) box = box.Extend(shapes[i].BoundingBox());
		return box;
	}
	static BoxForTriangles(shapes) {
		if (shapes.length === 0) return new ln.Box(new ln.Vector(0, 0, 0), new ln.Vector(0, 0, 0));
		let box = shapes[0].BoundingBox();
		for (let i = 1; i < shapes.length; i++) box.Extend(shapes[i].BoundingBox());
		return box;
	}
	static BoxForVectors(vectors) {
		if (vectors.length === 0) return new ln.Box(new ln.Vector(0, 0, 0), new ln.Vector(0, 0, 0));
		let [min, max] = [vectors[0], vectors[0]];
		for (let i = 1; i < vectors.length; i++) [min, max] = [min.Min(vectors[i]), max.Max(vectors[i])];
		return new ln.Box(min, max);
	}
	Size() { return this.Max.Sub(this.Min); }
	Anchor(anchor) { return this.Min.Add(this.Size().Mul(anchor)); }
	Center() { return this.Anchor(new ln.Vector(0.5, 0.5, 0.5)); }
	Contains(v) {
		return (
			this.Min.X <= v.X && this.Max.X >= v.X &&
			this.Min.Y <= v.Y && this.Max.Y >= v.Y &&
			this.Min.Z <= v.Z && this.Max.Z >= v.Z
		);
	}
	Extend(b) { return new ln.Box(this.Min.Min(b.Min), this.Max.Max(b.Max)); }
	Intersect(r) {
		let x1 = (this.Min.X - r.Origin.X) / r.Direction.X;
		let y1 = (this.Min.Y - r.Origin.Y) / r.Direction.Y;
		let z1 = (this.Min.Z - r.Origin.Z) / r.Direction.Z;
		let x2 = (this.Max.X - r.Origin.X) / r.Direction.X;
		let y2 = (this.Max.Y - r.Origin.Y) / r.Direction.Y;
		let z2 = (this.Max.Z - r.Origin.Z) / r.Direction.Z;
		if (x1 > x2) [x1, x2] = [x2, x1];
		if (y1 > y2) [y1, y2] = [y2, y1];
		if (z1 > z2) [z1, z2] = [z2, z1];
		return [Math.max(x1, y1, z1), Math.min(x2, y2, z2)]
	}
	Partition(axis, point) {
		if (axis === ln.AXIS.X) return [this.Min.X <= point, this.Max.X >= point];
		if (axis === ln.AXIS.Y) return [this.Min.Y <= point, this.Max.Y >= point];
		if (axis === ln.AXIS.Z) return [this.Min.Z <= point, this.Max.Z >= point];
		return [null, null];
	}
};

ln.Ray = class {
	constructor(origin, direction) { [this.Origin, this.Direction] = [origin, direction]; }
	Position(t) { return this.Origin.Add(this.Direction.MulScalar(t)); }
};

ln.Vector = class {
	constructor(x, y, z) { [this.X, this.Y, this.Z] = [x, y, z]; }
	static RandomUnitVector() {
		let x, y, z;
		while (true) {
			[x, y, z] = [2 * Math.Random() - 1, 2 * Math.Random() - 1, 2 * Math.Random() - 1];
			if (x * x + y * y + z * z > 1) break;
		}
		return new ln.Vector(x, y, z).Normalize();
	}
	Add(v) { return new ln.Vector(this.X + v.X, this.Y + v.Y, this.Z + v.Z); }
	Sub(v) { return new ln.Vector(this.X - v.X, this.Y - v.Y, this.Z - v.Z); }
	Mul(v) { return new ln.Vector(this.X * v.X, this.Y * v.Y, this.Z * v.Z); }
	Div(v) { return new ln.Vector(this.X / v.X, this.Y / v.Y, this.Z / v.Z); }
	AddScalar(s) { return new ln.Vector(this.X + s, this.Y + s, this.Z + s); }
	SubScalar(s) { return new ln.Vector(this.X - s, this.Y - s, this.Z - s); }
	MulScalar(s) { return new ln.Vector(this.X * s, this.Y * s, this.Z * s); }
	DivScalar(s) { return new ln.Vector(this.X / s, this.Y / s, this.Z / s); }
	Normalize() { return this.DivScalar(this.Length()); }
	LengthSquared() { return this.X * this.X + this.Y * this.Y + this.Z * this.Z; }
	Length() { return Math.sqrt(this.LengthSquared()); }
	DistanceSquared(v) { return this.Sub(v).LengthSquared(); }
	Distance(v) { return this.Sub(v).Length(); }
	Dot(v) { return this.X * v.X + this.Y * v.Y + this.Z * v.Z; }
	Cross(v) { return new ln.Vector(this.Y * v.Z - this.Z * v.Y, this.Z * v.X - this.X * v.Z, this.X * v.Y - this.Y * v.X); }
	Min(v) { return new ln.Vector(Math.min(this.X, v.X), Math.min(this.Y, v.Y), Math.min(this.Z, v.Z)); }
	Max(v) { return new ln.Vector(Math.max(this.X, v.X), Math.max(this.Y, v.Y), Math.max(this.Z, v.Z)); }
	MinComponent() { return Math.min(this.X, this.Y, this.Z); }
	MinAxis() {
		let [x, y, z] = [Math.abs(this.X), Math.abs(this.Y), Math.abs(this.Z)];
		if (x <= y && x <= z) return new ln.Vector(1, 0, 0);
		else if (y <= x && y <= z) return new ln.Vector(0, 1, 0);
		else return new ln.Vector(0, 0, 1);
	}
	SegmentDistance(v, w) {
		let l2 = v.DistanceSquared(w);
		if (l2 === 0) return this.Distance(v);
		let t = this.Sub(v).Dot(w.Sub(v)) / l2;
		if (t < 0) return this.Distance(v);
		else if (t > 1) return this.Distance(w);
		else return v.Add(w.Sub(v).MulScalar(t)).Distance(this);
	}
};

ln.Matrix = class {
	constructor(matrix) {
		[
			this.x00, this.x01, this.x02, this.x03,
			this.x10, this.x11, this.x12, this.x13,
			this.x20, this.x21, this.x22, this.x23,
			this.x30, this.x31, this.x32, this.x33
		] = matrix;
	}
	static createIdentity() {
		return new ln.Matrix([
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1
		]);
	}
	static createTranslate(v) {
		return new ln.Matrix([
			1, 0, 0, v.X,
			0, 1, 0, v.Y,
			0, 0, 1, v.Z,
			0, 0, 0,   1
		]);
	}
	static createScale(v) {
		return new ln.Matrix([
			v.X,   0,   0, 0,
			  0, v.Y,   0, 0,
			  0,   0, v.Z, 0,
			  0,   0,   0, 1
		]);
	}
	static createRotate(v, a) {
		v = v.Normalize();
		let [s, c] = [Math.sin(a), Math.cos(a)];
		let m = 1 - c;
		return new ln.Matrix([
			    m*v.X*v.X + c, m*v.X*v.Y + v.Z*s, m*v.Z*v.X - v.Y*s, 0,
			m*v.X*v.Y - v.Z*s,     m*v.Y*v.Y + c, m*v.Y*v.Z + v.X*s, 0,
			m*v.Z*v.X + v.Y*s, m*v.Y*v.Z - v.X*s,     m*v.Z*v.Z + c, 0,
			                0,                 0,                 0, 1
		]);
	}
	static createFrustum(l, r, b, t, n, f) {
		let [t1, t2, t3, t4] = [2 * n, r - l, t - b, f - n];
		return new ln.Matrix([
			t1 / t2,       0,  (r + l) / t2,              0,
			      0, t1 / t3,  (t + b) / t3,              0,
			      0,       0, (-f - n) / t4, (-t1 * f) / t4,
			      0,       0,            -1,              0
		]);
	}
	static createOrthographic(l, r, b, t, n, f) {
		return new ln.Matrix([
			2 / (r - l),           0,            0, -(r + l) / (r - l),
			          0, 2 / (t - b),            0, -(t + b) / (t - b),
			          0,           0, -2 / (f - n), -(f + n) / (f - n),
			          0,           0,            0,                  1
		]);
	}
	static createPerspective(fovy, aspect, near, far) {
		let ymax = near * Math.tan(fovy * Math.PI / 360);
		let xmax = ymax * aspect;
		return ln.Matrix.createFrustum(-xmax, xmax, -ymax, ymax, near, far);
	}
	static createLookAt(eye, center, up) {
		up = up.Normalize();
		let f = center.Sub(eye).Normalize();
		let s = f.Cross(up).Normalize();
		let u = s.Cross(f).Normalize();
		let m = new ln.Matrix([
			s.X, u.X, -f.X, eye.X,
			s.Y, u.Y, -f.Y, eye.Y,
			s.Z, u.Z, -f.Z, eye.Z,
			  0,   0,    0,     1
		]);
		return m.Inverse();
	}
	Translate(v) { return ln.Matrix.createTranslate(v).Mul(this); }
	Scale(v) { return ln.Matrix.createScale(v).Mul(this); }
	Rotate(v, a) { return ln.Matrix.createRotate(v, a).Mul(this); }
	Frustum(l, r, b, t, n, f) { return ln.Matrix.createFrustum(l, r, b, t, n, f).Mul(this); }
	Orthographic(l, r, b, t, n, f) { return ln.Matrix.createOrthographic(l, r, b, t, n, f).Mul(this); }
	Perspective(fovy, aspect, near, far) { return ln.Matrix.createPerspective(fovy, aspect, near, far).Mul(this); }
	Mul(m) {
		return new ln.Matrix([
			this.x00 * m.x00 + this.x01 * m.x10 + this.x02 * m.x20 + this.x03 * m.x30,
			this.x00 * m.x01 + this.x01 * m.x11 + this.x02 * m.x21 + this.x03 * m.x31,
			this.x00 * m.x02 + this.x01 * m.x12 + this.x02 * m.x22 + this.x03 * m.x32,
			this.x00 * m.x03 + this.x01 * m.x13 + this.x02 * m.x23 + this.x03 * m.x33,
			this.x10 * m.x00 + this.x11 * m.x10 + this.x12 * m.x20 + this.x13 * m.x30,
			this.x10 * m.x01 + this.x11 * m.x11 + this.x12 * m.x21 + this.x13 * m.x31,
			this.x10 * m.x02 + this.x11 * m.x12 + this.x12 * m.x22 + this.x13 * m.x32,
			this.x10 * m.x03 + this.x11 * m.x13 + this.x12 * m.x23 + this.x13 * m.x33,
			this.x20 * m.x00 + this.x21 * m.x10 + this.x22 * m.x20 + this.x23 * m.x30,
			this.x20 * m.x01 + this.x21 * m.x11 + this.x22 * m.x21 + this.x23 * m.x31,
			this.x20 * m.x02 + this.x21 * m.x12 + this.x22 * m.x22 + this.x23 * m.x32,
			this.x20 * m.x03 + this.x21 * m.x13 + this.x22 * m.x23 + this.x23 * m.x33,
			this.x30 * m.x00 + this.x31 * m.x10 + this.x32 * m.x20 + this.x33 * m.x30,
			this.x30 * m.x01 + this.x31 * m.x11 + this.x32 * m.x21 + this.x33 * m.x31,
			this.x30 * m.x02 + this.x31 * m.x12 + this.x32 * m.x22 + this.x33 * m.x32,
			this.x30 * m.x03 + this.x31 * m.x13 + this.x32 * m.x23 + this.x33 * m.x33
		]);
	}
	MulPosition(v) {
		let x = this.x00 * v.X + this.x01 * v.Y + this.x02 * v.Z + this.x03;
		let y = this.x10 * v.X + this.x11 * v.Y + this.x12 * v.Z + this.x13;
		let z = this.x20 * v.X + this.x21 * v.Y + this.x22 * v.Z + this.x23;
		return new ln.Vector(x, y, z);
	}
	MulPositionW(v) {
		let x = this.x00 * v.X + this.x01 * v.Y + this.x02 * v.Z + this.x03;
		let y = this.x10 * v.X + this.x11 * v.Y + this.x12 * v.Z + this.x13;
		let z = this.x20 * v.X + this.x21 * v.Y + this.x22 * v.Z + this.x23;
		let w = this.x30 * v.X + this.x31 * v.Y + this.x32 * v.Z + this.x33;
		return new ln.Vector(x / w, y / w, z / w);
	}
	MulDirection(v) {
		let x = this.x00 * v.X + this.x01 * v.Y + this.x02 * v.Z;
		let y = this.x10 * v.X + this.x11 * v.Y + this.x12 * v.Z;
		let z = this.x20 * v.X + this.x21 * v.Y + this.x22 * v.Z;
		return new ln.Vector(x, y, z).Normalize();
	}
	MulRay(r) {
		return new ln.Ray(this.MulPosition(r.Origin), this.MulDirection(r.Direction));
	}
	MulBox(box) {
		let r = new ln.Vector(this.x00, this.x10, this.x20);
		let u = new ln.Vector(this.x01, this.x11, this.x21);
		let b = new ln.Vector(this.x02, this.x12, this.x22);
		let t = new ln.Vector(this.x03, this.x13, this.x23);
		let xa = r.MulScalar(box.Min.X);
		let xb = r.MulScalar(box.Max.X);
		let ya = u.MulScalar(box.Min.Y);
		let yb = u.MulScalar(box.Max.Y);
		let za = b.MulScalar(box.Min.Z);
		let zb = b.MulScalar(box.Max.Z);
		[xa, xb] = [xa.Min(xb), xa.Max(xb)];
		[ya, yb] = [ya.Min(yb), ya.Max(yb)];
		[za, zb] = [za.Min(zb), za.Max(zb)];
		let min = xa.Add(ya).Add(za).Add(t);
		let max = xb.Add(yb).Add(zb).Add(t);
		return new ln.Box(min, max);
	}
	Transpose() {
		return new ln.Matrix([
			this.x00, this.x10, this.x20, this.x30,
			this.x01, this.x11, this.x21, this.x31,
			this.x02, this.x12, this.x22, this.x32,
			this.x03, this.x13, this.x23, this.x33
		]);
	}
	Determinant() {
		return (
			this.x00 * this.x11 * this.x22 * this.x33 - this.x00 * this.x11 * this.x23 * this.x32 +
			this.x00 * this.x12 * this.x23 * this.x31 - this.x00 * this.x12 * this.x21 * this.x33 +
			this.x00 * this.x13 * this.x21 * this.x32 - this.x00 * this.x13 * this.x22 * this.x31 -
			this.x01 * this.x12 * this.x23 * this.x30 + this.x01 * this.x12 * this.x20 * this.x33 -
			this.x01 * this.x13 * this.x20 * this.x32 + this.x01 * this.x13 * this.x22 * this.x30 -
			this.x01 * this.x10 * this.x22 * this.x33 + this.x01 * this.x10 * this.x23 * this.x32 +
			this.x02 * this.x13 * this.x20 * this.x31 - this.x02 * this.x13 * this.x21 * this.x30 +
			this.x02 * this.x10 * this.x21 * this.x33 - this.x02 * this.x10 * this.x23 * this.x31 +
			this.x02 * this.x11 * this.x23 * this.x30 - this.x02 * this.x11 * this.x20 * this.x33 -
			this.x03 * this.x10 * this.x21 * this.x32 + this.x03 * this.x10 * this.x22 * this.x31 -
			this.x03 * this.x11 * this.x22 * this.x30 + this.x03 * this.x11 * this.x20 * this.x32 -
			this.x03 * this.x12 * this.x20 * this.x31 + this.x03 * this.x12 * this.x21 * this.x30
		);
	}
	Inverse() {
		let d = this.Determinant();
		return new ln.Matrix([
			(this.x12 * this.x23 * this.x31 - this.x13 * this.x22 * this.x31 + this.x13 * this.x21 * this.x32 - this.x11 * this.x23 * this.x32 - this.x12 * this.x21 * this.x33 + this.x11 * this.x22 * this.x33) / d,
			(this.x03 * this.x22 * this.x31 - this.x02 * this.x23 * this.x31 - this.x03 * this.x21 * this.x32 + this.x01 * this.x23 * this.x32 + this.x02 * this.x21 * this.x33 - this.x01 * this.x22 * this.x33) / d,
			(this.x02 * this.x13 * this.x31 - this.x03 * this.x12 * this.x31 + this.x03 * this.x11 * this.x32 - this.x01 * this.x13 * this.x32 - this.x02 * this.x11 * this.x33 + this.x01 * this.x12 * this.x33) / d,
			(this.x03 * this.x12 * this.x21 - this.x02 * this.x13 * this.x21 - this.x03 * this.x11 * this.x22 + this.x01 * this.x13 * this.x22 + this.x02 * this.x11 * this.x23 - this.x01 * this.x12 * this.x23) / d,
			(this.x13 * this.x22 * this.x30 - this.x12 * this.x23 * this.x30 - this.x13 * this.x20 * this.x32 + this.x10 * this.x23 * this.x32 + this.x12 * this.x20 * this.x33 - this.x10 * this.x22 * this.x33) / d,
			(this.x02 * this.x23 * this.x30 - this.x03 * this.x22 * this.x30 + this.x03 * this.x20 * this.x32 - this.x00 * this.x23 * this.x32 - this.x02 * this.x20 * this.x33 + this.x00 * this.x22 * this.x33) / d,
			(this.x03 * this.x12 * this.x30 - this.x02 * this.x13 * this.x30 - this.x03 * this.x10 * this.x32 + this.x00 * this.x13 * this.x32 + this.x02 * this.x10 * this.x33 - this.x00 * this.x12 * this.x33) / d,
			(this.x02 * this.x13 * this.x20 - this.x03 * this.x12 * this.x20 + this.x03 * this.x10 * this.x22 - this.x00 * this.x13 * this.x22 - this.x02 * this.x10 * this.x23 + this.x00 * this.x12 * this.x23) / d,
			(this.x11 * this.x23 * this.x30 - this.x13 * this.x21 * this.x30 + this.x13 * this.x20 * this.x31 - this.x10 * this.x23 * this.x31 - this.x11 * this.x20 * this.x33 + this.x10 * this.x21 * this.x33) / d,
			(this.x03 * this.x21 * this.x30 - this.x01 * this.x23 * this.x30 - this.x03 * this.x20 * this.x31 + this.x00 * this.x23 * this.x31 + this.x01 * this.x20 * this.x33 - this.x00 * this.x21 * this.x33) / d,
			(this.x01 * this.x13 * this.x30 - this.x03 * this.x11 * this.x30 + this.x03 * this.x10 * this.x31 - this.x00 * this.x13 * this.x31 - this.x01 * this.x10 * this.x33 + this.x00 * this.x11 * this.x33) / d,
			(this.x03 * this.x11 * this.x20 - this.x01 * this.x13 * this.x20 - this.x03 * this.x10 * this.x21 + this.x00 * this.x13 * this.x21 + this.x01 * this.x10 * this.x23 - this.x00 * this.x11 * this.x23) / d,
			(this.x12 * this.x21 * this.x30 - this.x11 * this.x22 * this.x30 - this.x12 * this.x20 * this.x31 + this.x10 * this.x22 * this.x31 + this.x11 * this.x20 * this.x32 - this.x10 * this.x21 * this.x32) / d,
			(this.x01 * this.x22 * this.x30 - this.x02 * this.x21 * this.x30 + this.x02 * this.x20 * this.x31 - this.x00 * this.x22 * this.x31 - this.x01 * this.x20 * this.x32 + this.x00 * this.x21 * this.x32) / d,
			(this.x02 * this.x11 * this.x30 - this.x01 * this.x12 * this.x30 - this.x02 * this.x10 * this.x31 + this.x00 * this.x12 * this.x31 + this.x01 * this.x10 * this.x32 - this.x00 * this.x11 * this.x32) / d,
			(this.x01 * this.x12 * this.x20 - this.x02 * this.x11 * this.x20 + this.x02 * this.x10 * this.x21 - this.x00 * this.x12 * this.x21 - this.x01 * this.x10 * this.x22 + this.x00 * this.x11 * this.x22) / d,
		]);
	}
};

ln.Function = class {
	constructor(func, box, direction) { [this.Func, this.Box, this.Direction] = [func, box, direction]; }
	static Below() { return "Below"; }
	static Above() { return "Above"; }
	Compile() { return null; }
	BoundingBox() { return this.Box; }
	Contains(v) {
		if (this.Direction === "Below") return v.Z < this.Func(v.X, v.Y);
		else return v.Z > this.Func(v.X, v.Y);
	}
	Intersect(ray) {
		let step = 1 / 64;
		let sign = this.Contains(ray.Position(step));
		for (let t = step; t < 10; t += step) {
			let v = ray.Position(t);
			if (this.Contains(v) !== sign && this.Box.Contains(v)) return new ln.Hit(this, t);
		}
		return ln.Hit.none();
	}
	Paths() {
		let paths = new ln.Paths();
		let fine = 1 / 256;
		for (let a = 0; a < 360; a += 5) {
			let path = new ln.Path();
			for (let r = 0; r <= 8; r += fine) {
				let [x, y] = [r * Math.cos(ln.Utils.degToRad(a)), r * Math.sin(ln.Utils.degToRad(a))];
				let z = this.Func(x, y);
				let o = -Math.pow(Math.abs(-z), 1.4);
				o = 0;
				[x, y] = [r * (Math.cos(ln.Utils.degToRad(a) - o)), r * (Math.sin(ln.Utils.degToRad(a) - o))];
				z = Math.min(z, this.Box.Max.Z);
				z = Math.max(z, this.Box.Min.Z);
				path.Add(new ln.Vector(x, y, z));
			}
			paths.Add(path);
		}

		for (let r = 0.1; r < 10; r += 0.15) {
			let path = new ln.Path();
			for (let a = 0; a <= 360; a += 5) {
				let [x, y] = [r * Math.cos(ln.Utils.degToRad(a)), r * Math.sin(ln.Utils.degToRad(a))];
				let z = this.Func(x, y);
				[x, y] = [r * Math.cos(ln.Utils.degToRad(a)), r * Math.sin(ln.Utils.degToRad(a))];
				z = Math.min(z, this.Box.Max.Z);
				z = Math.max(z, this.Box.Min.Z);
				path.Add(new ln.Vector(x, y, z));
			}
			paths.Add(path);
		}

		return paths;
	}
	// Paths() {
	// 	let paths = new ln.Paths();
	// 	let [step, fine] = [1 / 16, 1 / 64];
	// 	for (let x = this.Box.Min.X; x <= this.Box.Max.X; x += step) {
	// 		let path = new ln.Path();
	// 		for (let y = this.Box.Min.Y; y <= this.Box.Max.Y; y += fine) {
	// 			let z = this.Func(x, y);
	// 			z = Math.min(z, this.Box.Max.Z);
	// 			z = Math.max(z, this.Box.Min.Z);
	// 			path.Add(new ln.Vector(x, y, z));
	// 		}
	// 		paths.Add(path);
	// 	}
	// 	for (let y = this.Box.Min.Y; y <= this.Box.Max.Y; y += step) {
	// 		let path = new ln.Path();
	// 		for (let x = this.Box.Min.X; x <= this.Box.Max.X; x += fine) {
	// 			let z = this.Func(x, y);
	// 			z = Math.min(z, this.Box.Max.Z);
	// 			z = Math.max(z, this.Box.Min.Z);
	// 			path.Add(new ln.Vector(x, y, z));
	// 		}
	// 		paths.Add(path);
	// 	}
	// 	return paths;
	// }
};

ln.Cube = class {
	constructor(min, max) {
		[this.Min, this.Max] = [min, max];
		this.Box = new ln.Box(this.Min, this.Max);
	}
	Compile() { return; }
	BoundingBox() { return this.Box; }
	Contains(v, f) {
		if (v.X < this.Min.X - f || v.X > this.Max.X + f) return false;
		if (v.Y < this.Min.Y - f || v.Y > this.Max.Y + f) return false;
		if (v.Z < this.Min.Z - f || v.Z > this.Max.Z + f) return false;
		return true;
	}
	Intersect(r) {
		let [n, f] = [this.Min.Sub(r.Origin).Div(r.Direction), this.Max.Sub(r.Origin).Div(r.Direction)];
		[n, f] = [n.Min(f), n.Max(f)];
		let [t0, t1] = [Math.max(n.X, n.Y, n.Z), Math.min(f.X, f.Y, f.Z)];
		if (t0 < 0.001 && t1 > 0.001) return new ln.Hit(this, t1);
		if (t0 >= 0.001 && t0 < t1) return new ln.Hit(this, t0);
		return ln.Hit.none();
	}
	Paths() {
		let [x1, y1, z1] = [this.Min.X, this.Min.Y, this.Min.Z];
		let [x2, y2, z2] = [this.Max.X, this.Max.Y, this.Max.Z];
		let paths = new ln.Paths();
		paths.Add(ln.Path.FromPoints([new ln.Vector(x1, y1, z1), new ln.Vector(x1, y1, z2)]));
		paths.Add(ln.Path.FromPoints([new ln.Vector(x1, y1, z1), new ln.Vector(x1, y2, z1)]));
		paths.Add(ln.Path.FromPoints([new ln.Vector(x1, y1, z1), new ln.Vector(x2, y1, z1)]));
		paths.Add(ln.Path.FromPoints([new ln.Vector(x1, y1, z2), new ln.Vector(x1, y2, z2)]));
		paths.Add(ln.Path.FromPoints([new ln.Vector(x1, y1, z2), new ln.Vector(x2, y1, z2)]));
		paths.Add(ln.Path.FromPoints([new ln.Vector(x1, y2, z1), new ln.Vector(x1, y2, z2)]));
		paths.Add(ln.Path.FromPoints([new ln.Vector(x1, y2, z1), new ln.Vector(x2, y2, z1)]));
		paths.Add(ln.Path.FromPoints([new ln.Vector(x1, y2, z2), new ln.Vector(x2, y2, z2)]));
		paths.Add(ln.Path.FromPoints([new ln.Vector(x2, y1, z1), new ln.Vector(x2, y1, z2)]));
		paths.Add(ln.Path.FromPoints([new ln.Vector(x2, y1, z1), new ln.Vector(x2, y2, z1)]));
		paths.Add(ln.Path.FromPoints([new ln.Vector(x2, y1, z2), new ln.Vector(x2, y2, z2)]));
		paths.Add(ln.Path.FromPoints([new ln.Vector(x2, y2, z1), new ln.Vector(x2, y2, z2)]));
		return paths;
	}
};







