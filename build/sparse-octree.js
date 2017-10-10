/**
 * sparse-octree v4.0.2 build Oct 10 2017
 * https://github.com/vanruesc/sparse-octree
 * Copyright 2017 Raoul van Rüschen, Zlib
 */

(function (global, factory) {
	typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
	typeof define === 'function' && define.amd ? define(['exports'], factory) :
	(factory((global.SPARSEOCTREE = {})));
}(this, (function (exports) { 'use strict';

	/**
	 * A vector with three components.
	 */

	var Vector3 = function Vector3(x, y, z) {
		if ( x === void 0 ) x = 0;
		if ( y === void 0 ) y = 0;
		if ( z === void 0 ) z = 0;


		/**
			 * The X component.
			 *
			 * @type {Number}
			 */

		this.x = x;

		/**
			 * The Y component.
			 *
			 * @type {Number}
			 */

		this.y = y;

		/**
			 * The Z component.
			 *
			 * @type {Number}
			 */

		this.z = z;

	};

	/**
		 * Sets the values of this vector
		 *
		 * @param {Number} x - The X component.
		 * @param {Number} y - The Y component.
		 * @param {Number} z - The Z component.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.set = function set (x, y, z) {

		this.x = x;
		this.y = y;
		this.z = z;

		return this;

	};

	/**
		 * Copies the values of another vector.
		 *
		 * @param {Vector3} v - A vector.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.copy = function copy (v) {

		this.x = v.x;
		this.y = v.y;
		this.z = v.z;

		return this;

	};

	/**
		 * Clones this vector.
		 *
		 * @return {Vector3} A clone of this vector.
		 */

	Vector3.prototype.clone = function clone () {

		return new this.constructor(this.x, this.y, this.z);

	};

	/**
		 * Copies values from an array.
		 *
		 * @param {Number[]} array - An array.
		 * @param {Number} offset - An offset.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.fromArray = function fromArray (array, offset) {
			if ( offset === void 0 ) offset = 0;


		this.x = array[offset];
		this.y = array[offset + 1];
		this.z = array[offset + 2];

		return this;

	};

	/**
		 * Stores this vector in an array.
		 *
		 * @param {Array} [array] - A target array.
		 * @param {Number} offset - An offset.
		 * @return {Number[]} The array.
		 */

	Vector3.prototype.toArray = function toArray (array, offset) {
			if ( array === void 0 ) array = [];
			if ( offset === void 0 ) offset = 0;


		array[offset] = this.x;
		array[offset + 1] = this.y;
		array[offset + 2] = this.z;

		return array;

	};

	/**
		 * Sets the values of this vector based on a spherical description.
		 *
		 * @param {Spherical} s - A spherical description.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.setFromSpherical = function setFromSpherical (s) {

		var sinPhiRadius = Math.sin(s.phi) * s.radius;

		this.x = sinPhiRadius * Math.sin(s.theta);
		this.y = Math.cos(s.phi) * s.radius;
		this.z = sinPhiRadius * Math.cos(s.theta);

		return this;

	};

	/**
		 * Sets the values of this vector based on a cylindrical description.
		 *
		 * @param {Cylindrical} c - A cylindrical description.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.setFromCylindrical = function setFromCylindrical (c) {

		this.x = c.radius * Math.sin(c.theta);
		this.y = c.y;
		this.z = c.radius * Math.cos(c.theta);

		return this;

	};

	/**
		 * Copies the values of a matrix column.
		 *
		 * @param {Matrix4} m - A 4x4 matrix.
		 * @param {Number} index - A column index of the range [0, 2].
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.setFromMatrixColumn = function setFromMatrixColumn (m, index) {

		return this.fromArray(m.elements, index * 4);

	};

	/**
		 * Extracts the position from a matrix.
		 *
		 * @param {Matrix4} m - A 4x4 matrix.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.setFromMatrixPosition = function setFromMatrixPosition (m) {

		var me = m.elements;

		this.x = me[12];
		this.y = me[13];
		this.z = me[14];

		return this;

	};

	/**
		 * Extracts the scale from a matrix.
		 *
		 * @param {Matrix4} m - A 4x4 matrix.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.setFromMatrixScale = function setFromMatrixScale (m) {

		var sx = this.setFromMatrixColumn(m, 0).length();
		var sy = this.setFromMatrixColumn(m, 1).length();
		var sz = this.setFromMatrixColumn(m, 2).length();

		this.x = sx;
		this.y = sy;
		this.z = sz;

		return this;

	};

	/**
		 * Adds a vector to this one.
		 *
		 * @param {Vector3} v - The vector to add.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.add = function add (v) {

		this.x += v.x;
		this.y += v.y;
		this.z += v.z;

		return this;

	};

	/**
		 * Adds a scalar to this vector.
		 *
		 * @param {Number} s - The scalar to add.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.addScalar = function addScalar (s) {

		this.x += s;
		this.y += s;
		this.z += s;

		return this;

	};

	/**
		 * Sets this vector to the sum of two given vectors.
		 *
		 * @param {Vector3} a - A vector.
		 * @param {Vector3} b - Another vector.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.addVectors = function addVectors (a, b) {

		this.x = a.x + b.x;
		this.y = a.y + b.y;
		this.z = a.z + b.z;

		return this;

	};

	/**
		 * Adds a scaled vector to this one.
		 *
		 * @param {Vector3} v - The vector to scale and add.
		 * @param {Number} s - A scalar.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.addScaledVector = function addScaledVector (v, s) {

		this.x += v.x * s;
		this.y += v.y * s;
		this.z += v.z * s;

		return this;

	};

	/**
		 * Subtracts a vector from this vector.
		 *
		 * @param {Vector3} v - The vector to subtract.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.sub = function sub (v) {

		this.x -= v.x;
		this.y -= v.y;
		this.z -= v.z;

		return this;

	};

	/**
		 * Subtracts a scalar to this vector.
		 *
		 * @param {Number} s - The scalar to subtract.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.subScalar = function subScalar (s) {

		this.x -= s;
		this.y -= s;
		this.z -= s;

		return this;

	};

	/**
		 * Sets this vector to the difference between two given vectors.
		 *
		 * @param {Vector3} a - A vector.
		 * @param {Vector3} b - A second vector.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.subVectors = function subVectors (a, b) {

		this.x = a.x - b.x;
		this.y = a.y - b.y;
		this.z = a.z - b.z;

		return this;

	};

	/**
		 * Multiplies this vector with another vector.
		 *
		 * @param {Vector3} v - A vector.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.multiply = function multiply (v) {

		this.x *= v.x;
		this.y *= v.y;
		this.z *= v.z;

		return this;

	};

	/**
		 * Multiplies this vector with a given scalar.
		 *
		 * @param {Number} s - A scalar.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.multiplyScalar = function multiplyScalar (s) {

		this.x *= s;
		this.y *= s;
		this.z *= s;

		return this;

	};

	/**
		 * Sets this vector to the product of two given vectors.
		 *
		 * @param {Vector3} a - A vector.
		 * @param {Vector3} b - Another vector.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.multiplyVectors = function multiplyVectors (a, b) {

		this.x = a.x * b.x;
		this.y = a.y * b.y;
		this.z = a.z * b.z;

		return this;

	};

	/**
		 * Divides this vector by another vector.
		 *
		 * @param {Vector3} v - A vector.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.divide = function divide (v) {

		this.x /= v.x;
		this.y /= v.y;
		this.z /= v.z;

		return this;

	};

	/**
		 * Divides this vector by a given scalar.
		 *
		 * @param {Number} s - A scalar.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.divideScalar = function divideScalar (s) {

		this.x /= s;
		this.y /= s;
		this.z /= s;

		return this;

	};

	/**
		 * Calculates the cross product of this vector and the given one.
		 *
		 * @param {Vector3} v - A vector.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.cross = function cross (v) {

		var x = this.x, y = this.y, z = this.z;

		this.x = y * v.z - z * v.y;
		this.y = z * v.x - x * v.z;
		this.z = x * v.y - y * v.x;

		return this;

	};

	/**
		 * Sets this vector to the cross product of the given vectors.
		 *
		 * @param {Vector3} a - A vector.
		 * @param {Vector3} b - Another vector.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.crossVectors = function crossVectors (a, b) {

		var ax = a.x, ay = a.y, az = a.z;
		var bx = b.x, by = b.y, bz = b.z;

		this.x = ay * bz - az * by;
		this.y = az * bx - ax * bz;
		this.z = ax * by - ay * bx;

		return this;

	};

	/**
		 * Applies a matrix to this vector.
		 *
		 * @param {Matrix3} m - A matrix.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.applyMatrix3 = function applyMatrix3 (m) {

		var x = this.x, y = this.y, z = this.z;
		var e = m.elements;

		this.x = e[0] * x + e[3] * y + e[6] * z;
		this.y = e[1] * x + e[4] * y + e[7] * z;
		this.z = e[2] * x + e[5] * y + e[8] * z;

		return this;

	};

	/**
		 * Applies a matrix to this vector.
		 *
		 * @param {Matrix4} m - A matrix.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.applyMatrix4 = function applyMatrix4 (m) {

		var x = this.x, y = this.y, z = this.z;
		var e = m.elements;

		this.x = e[0] * x + e[4] * y + e[8] * z + e[12];
		this.y = e[1] * x + e[5] * y + e[9] * z + e[13];
		this.z = e[2] * x + e[6] * y + e[10] * z + e[14];

		return this;

	};

	/**
		 * Applies a quaternion to this vector.
		 *
		 * @param {Quaternion} q - A quaternion.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.applyQuaternion = function applyQuaternion (q) {

		var x = this.x, y = this.y, z = this.z;
		var qx = q.x, qy = q.y, qz = q.z, qw = q.w;

		// Calculate: quaternion * vector.
		var ix = qw * x + qy * z - qz * y;
		var iy = qw * y + qz * x - qx * z;
		var iz = qw * z + qx * y - qy * x;
		var iw = -qx * x - qy * y - qz * z;

		// Calculate: result * inverse quaternion.
		this.x = ix * qw + iw * -qx + iy * -qz - iz * -qy;
		this.y = iy * qw + iw * -qy + iz * -qx - ix * -qz;
		this.z = iz * qw + iw * -qz + ix * -qy - iy * -qx;

		return this;

	};

	/**
		 * Negates this vector.
		 *
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.negate = function negate () {

		this.x = -this.x;
		this.y = -this.y;
		this.z = -this.z;

		return this;

	};

	/**
		 * Calculates the dot product with another vector.
		 *
		 * @param {Vector3} v - A vector.
		 * @return {Number} The dot product.
		 */

	Vector3.prototype.dot = function dot (v) {

		return this.x * v.x + this.y * v.y + this.z * v.z;

	};

	/**
		 * Reflects this vector. The given plane normal is assumed to be normalized.
		 *
		 * @param {Vector3} n - A normal.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.reflect = function reflect (n, target) {
			if ( target === void 0 ) target = new Vector3();


		var nx = n.x;
		var ny = n.y;
		var nz = n.z;

		this.sub(n.multiplyScalar(2 * this.dot(n)));

		// Restore the normal.
		n.set(nx, ny, nz);

		return this;

	};

	/**
		 * Computes the angle to the given vector.
		 *
		 * @param {Vector3} v - A vector.
		 * @return {Number} The angle in radians.
		 */

	Vector3.prototype.angleTo = function angleTo (v) {

		var theta = this.dot(v) / (Math.sqrt(this.lengthSquared() * v.lengthSquared()));

		// Clamp to avoid numerical problems.
		return Math.acos(Math.min(Math.max(theta, -1), 1));

	};

	/**
		 * Calculates the Manhattan length of this vector.
		 *
		 * @return {Number} The length.
		 */

	Vector3.prototype.lengthManhattan = function lengthManhattan () {

		return Math.abs(this.x) + Math.abs(this.y) + Math.abs(this.z);

	};

	/**
		 * Calculates the squared length of this vector.
		 *
		 * @return {Number} The squared length.
		 */

	Vector3.prototype.lengthSquared = function lengthSquared () {

		return this.x * this.x + this.y * this.y + this.z * this.z;

	};

	/**
		 * Calculates the length of this vector.
		 *
		 * @return {Number} The length.
		 */

	Vector3.prototype.length = function length () {

		return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);

	};

	/**
		 * Calculates the Manhattan distance to a given vector.
		 *
		 * @param {Vector3} v - A vector.
		 * @return {Number} The distance.
		 */

	Vector3.prototype.distanceToManhattan = function distanceToManhattan (v) {

		return Math.abs(this.x - v.x) + Math.abs(this.y - v.y) + Math.abs(this.z - v.z);

	};

	/**
		 * Calculates the squared distance to a given vector.
		 *
		 * @param {Vector3} v - A vector.
		 * @return {Number} The squared distance.
		 */

	Vector3.prototype.distanceToSquared = function distanceToSquared (v) {

		var dx = this.x - v.x;
		var dy = this.y - v.y;
		var dz = this.z - v.z;

		return dx * dx + dy * dy + dz * dz;

	};

	/**
		 * Calculates the distance to a given vector.
		 *
		 * @param {Vector3} v - A vector.
		 * @return {Number} The distance.
		 */

	Vector3.prototype.distanceTo = function distanceTo (v) {

		return Math.sqrt(this.distanceToSquared(v));

	};

	/**
		 * Normalizes this vector.
		 *
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.normalize = function normalize () {

		return this.divideScalar(this.length());

	};

	/**
		 * Sets the length of this vector.
		 *
		 * @param {Number} length - The new length.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.setLength = function setLength (length) {

		return this.normalize().multiplyScalar(length);

	};

	/**
		 * Adopts the min value for each component of this vector and the given one.
		 *
		 * @param {Vector3} v - A vector.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.min = function min (v) {

		this.x = Math.min(this.x, v.x);
		this.y = Math.min(this.y, v.y);
		this.z = Math.min(this.z, v.z);

		return this;

	};

	/**
		 * Adopts the max value for each component of this vector and the given one.
		 *
		 * @param {Vector3} v - A vector.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.max = function max (v) {

		this.x = Math.max(this.x, v.x);
		this.y = Math.max(this.y, v.y);
		this.z = Math.max(this.z, v.z);

		return this;

	};

	/**
		 * Clamps this vector.
		 *
		 * @param {Vector3} min - The lower bounds. Assumed to be smaller than max.
		 * @param {Vector3} max - The upper bounds. Assumed to be greater than min.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.clamp = function clamp (min, max) {

		this.x = Math.max(min.x, Math.min(max.x, this.x));
		this.y = Math.max(min.y, Math.min(max.y, this.y));
		this.z = Math.max(min.z, Math.min(max.z, this.z));

		return this;

	};

	/**
		 * Floors this vector.
		 *
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.floor = function floor () {

		this.x = Math.floor(this.x);
		this.y = Math.floor(this.y);
		this.z = Math.floor(this.z);

		return this;

	};

	/**
		 * Ceils this vector.
		 *
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.ceil = function ceil () {

		this.x = Math.ceil(this.x);
		this.y = Math.ceil(this.y);
		this.z = Math.ceil(this.z);

		return this;

	};

	/**
		 * Rounds this vector.
		 *
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.round = function round () {

		this.x = Math.round(this.x);
		this.y = Math.round(this.y);
		this.z = Math.round(this.z);

		return this;

	};

	/**
		 * Lerps towards the given vector.
		 *
		 * @param {Vector3} v - The target vector.
		 * @param {Number} alpha - The lerp factor.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.lerp = function lerp (v, alpha) {

		this.x += (v.x - this.x) * alpha;
		this.y += (v.y - this.y) * alpha;
		this.z += (v.z - this.z) * alpha;

		return this;

	};

	/**
		 * Sets this vector to the lerp result of the given vectors.
		 *
		 * @param {Vector3} v1 - A base vector.
		 * @param {Vector3} v2 - The target vector.
		 * @param {Number} alpha - The lerp factor.
		 * @return {Vector3} This vector.
		 */

	Vector3.prototype.lerpVectors = function lerpVectors (v1, v2, alpha) {

		return this.subVectors(v2, v1).multiplyScalar(alpha).add(v1);

	};

	/**
		 * Checks if this vector equals the given one.
		 *
		 * @param {Vector3} v - A vector.
		 * @return {Boolean} Whether this vector equals the given one.
		 */

	Vector3.prototype.equals = function equals (v) {

		return (v.x === this.x && v.y === this.y && v.z === this.z);

	};

	/**
	 * A vector.
	 *
	 * @type {Vector3}
	 * @private
	 */

	var v$1 = new Vector3();

	/**
	 * A 3D box.
	 */

	var Box3 = function Box3(
		min,
		max
	) {
		if ( min === void 0 ) min = new Vector3(Infinity, Infinity, Infinity);
		if ( max === void 0 ) max = new Vector3(-Infinity, -Infinity, -Infinity);


		/**
			 * The lower bounds.
			 *
			 * @type {Vector3}
			 */

		this.min = min;

		/**
			 * The upper bounds.
			 *
			 * @type {Vector3}
			 */

		this.max = max;

	};

	/**
		 * Sets the values of this box.
		 *
		 * @param {Vector3} min - The lower bounds.
		 * @param {Vector3} max - The upper bounds.
		 * @return {Box3} This box.
		 */

	Box3.prototype.set = function set (min, max) {

		this.min.copy(min);
		this.max.copy(max);

		return this;

	};

	/**
		 * Copies the values of a given box.
		 *
		 * @param {Box3} b - A box.
		 * @return {Box3} This box.
		 */

	Box3.prototype.copy = function copy (b) {

		this.min.copy(b.min);
		this.max.copy(b.max);

		return this;

	};

	/**
		 * Clones this box.
		 *
		 * @return {Box3} A clone of this box.
		 */

	Box3.prototype.clone = function clone () {

		return new this.constructor().copy(this);

	};

	/**
		 * Makes this box empty.
		 *
		 * The lower bounds are set to infinity and the upper bounds to negative
		 * infinity to create an infinitely small box.
		 *
		 * @return {Box3} This box.
		 */

	Box3.prototype.makeEmpty = function makeEmpty () {

		this.min.x = this.min.y = this.min.z = Infinity;
		this.max.x = this.max.y = this.max.z = -Infinity;

		return this;

	};

	/**
		 * Indicates whether this box is truly empty.
		 *
		 * This is a more robust check for emptiness since the volume can get positive
		 * with two negative axes.
		 *
		 * @return {Box3} This box.
		 */

	Box3.prototype.isEmpty = function isEmpty () {

		return (
			this.max.x < this.min.x ||
			this.max.y < this.min.y ||
			this.max.z < this.min.z
		);

	};

	/**
		 * Computes the center of this box.
		 *
		 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
		 * @return {Vector3} A vector that describes the center of this box.
		 */

	Box3.prototype.getCenter = function getCenter (target) {
			if ( target === void 0 ) target = new Vector3();


		return !this.isEmpty() ?
			target.addVectors(this.min, this.max).multiplyScalar(0.5) :
			target.set(0, 0, 0);

	};

	/**
		 * Computes the size of this box.
		 *
		 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
		 * @return {Vector3} A vector that describes the size of this box.
		 */

	Box3.prototype.getSize = function getSize (target) {
			if ( target === void 0 ) target = new Vector3();


		return !this.isEmpty() ?
			target.subVectors(this.max, this.min) :
			target.set(0, 0, 0);

	};

	/**
		 * Computes the bounding sphere of this box.
		 *
		 * @param {Sphere} [target] - A target sphere. If none is provided, a new one will be created.
		 * @return {Sphere} The bounding sphere of this box.
		 */

	Box3.prototype.getBoundingSphere = function getBoundingSphere (target) {
			if ( target === void 0 ) target = new Sphere();


		this.getCenter(target.center);

		target.radius = this.getSize(v$1).length() * 0.5;

		return target;

	};

	/**
		 * Expands this box by the given point.
		 *
		 * @param {Vector3} p - A point.
		 * @return {Box3} This box.
		 */

	Box3.prototype.expandByPoint = function expandByPoint (p) {

		this.min.min(p);
		this.max.max(p);

		return this;

	};

	/**
		 * Expands this box by the given vector.
		 *
		 * @param {Vector3} v - A vector.
		 * @return {Box3} This box.
		 */

	Box3.prototype.expandByVector = function expandByVector (v) {

		this.min.sub(v);
		this.max.add(v);

		return this;

	};

	/**
		 * Expands this box by the given scalar.
		 *
		 * @param {Number} s - A scalar.
		 * @return {Box3} This box.
		 */

	Box3.prototype.expandByScalar = function expandByScalar (s) {

		this.min.addScalar(-s);
		this.max.addScalar(s);

		return this;

	};

	/**
		 * Defines this box by the given points.
		 *
		 * @param {Vector3[]} points - The points.
		 * @return {Box3} This box.
		 */

	Box3.prototype.setFromPoints = function setFromPoints (points) {
			var this$1 = this;


		var i, l;

		this.min.set(0, 0, 0);
		this.max.set(0, 0, 0);

		for(i = 0, l = points.length; i < l; ++i) {

			this$1.expandByPoint(points[i]);

		}

		return this;

	};

	/**
		 * Defines this box by the given center and size.
		 *
		 * @param {Vector3} center - The center.
		 * @param {Number} size - The size.
		 * @return {Box3} This box.
		 */

	Box3.prototype.setFromCenterAndSize = function setFromCenterAndSize (center, size) {

		var halfSize = v$1.copy(size).multiplyScalar(0.5);

		this.min.copy(center).sub(halfSize);
		this.max.copy(center).add(halfSize);

		return this;

	};

	/**
		 * Clamps the given point to the boundaries of this box.
		 *
		 * @param {Vector3} p - A point.
		 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
		 * @return {Vector3} The clamped point.
		 */

	Box3.prototype.clampPoint = function clampPoint (point, target) {
			if ( target === void 0 ) target = new Vector3();


		return target.copy(point).clamp(this.min, this.max);

	};

	/**
		 * Calculates the distance from this box to the given point.
		 *
		 * @param {Vector3} p - A point.
		 * @return {Number} The distance.
		 */

	Box3.prototype.distanceToPoint = function distanceToPoint (p) {

		var clampedPoint = v$1.copy(p).clamp(this.min, this.max);

		return clampedPoint.sub(p).length();

	};

	/**
		 * Translates this box.
		 *
		 * @param {Vector3} offset - The offset.
		 * @return {Box3} This box.
		 */

	Box3.prototype.translate = function translate (offset) {

		this.min.add(offset);
		this.max.add(offset);

		return this;

	};

	/**
		 * Expands this box by combining it with the given one.
		 *
		 * @param {Box3} b - A box.
		 * @return {Box3} This box.
		 */

	Box3.prototype.intersect = function intersect (b) {

		this.min.max(b.min);
		this.max.min(b.max);

		/* Ensure that if there is no overlap, the result is fully empty to prevent
		subsequent intersections to erroneously return valid values. */
		if(this.isEmpty()) { this.makeEmpty(); }

		return this;

	};

	/**
		 * Expands this box by combining it with the given one.
		 *
		 * @param {Box3} b - A box.
		 * @return {Box3} This box.
		 */

	Box3.prototype.union = function union (b) {

		this.min.min(b.min);
		this.max.max(b.max);

		return this;

	};

	/**
		 * Checks if the given point lies inside this box.
		 *
		 * @param {Vector3} p - A point.
		 * @return {Boolean} Whether this box contains the point.
		 */

	Box3.prototype.containsPoint = function containsPoint (p) {

		return !(
			p.x < this.min.x || p.x > this.max.x ||
			p.y < this.min.y || p.y > this.max.y ||
			p.z < this.min.z || p.z > this.max.z
		);

	};

	/**
		 * Checks if the given box lies inside this box.
		 *
		 * @param {Vector3} b - A box.
		 * @return {Boolean} Whether this box contains the given one.
		 */

	Box3.prototype.containsBox = function containsBox (b) {

		return (
			this.min.x <= b.min.x && b.max.x <= this.max.x &&
			this.min.y <= b.min.y && b.max.y <= this.max.y &&
			this.min.z <= b.min.z && b.max.z <= this.max.z
		);

	};

	/**
		 * Checks if this box intersects with the given one.
		 *
		 * @param {Box3} b - A box.
		 * @return {Boolean} Whether the boxes intersect.
		 */

	Box3.prototype.intersectsBox = function intersectsBox (b) {

		return !(
			b.max.x < this.min.x || b.min.x > this.max.x ||
			b.max.y < this.min.y || b.min.y > this.max.y ||
			b.max.z < this.min.z || b.min.z > this.max.z
		);

	};

	/**
		 * Checks if this box intersects with the given sphere.
		 *
		 * @param {Sphere} s - A sphere.
		 * @return {Boolean} Whether the box intersects with the sphere.
		 */

	Box3.prototype.intersectsSphere = function intersectsSphere (s) {

		// Find the point in this box that is closest to the sphere's center.
		var closestPoint = this.clampPoint(s.center, v$1);

		// If that point is inside the sphere, it intersects with this box.
		return (closestPoint.distanceToSquared(s.center) <= (s.radius * s.radius));

	};

	/**
		 * Checks if this box intersects with the given plane.
		 *
		 * Computes the minimum and maximum dot product values. If those values are on
		 * the same side (back or front) of the plane, then there is no intersection.
		 *
		 * @param {Plane} p - A plane.
		 * @return {Boolean} Whether the box intersects with the plane.
		 */

	Box3.prototype.intersectsPlane = function intersectsPlane (p) {

		var min, max;

		if(p.normal.x > 0) {

			min = p.normal.x * this.min.x;
			max = p.normal.x * this.max.x;

		} else {

			min = p.normal.x * this.max.x;
			max = p.normal.x * this.min.x;

		}

		if(p.normal.y > 0) {

			min += p.normal.y * this.min.y;
			max += p.normal.y * this.max.y;

		} else {

			min += p.normal.y * this.max.y;
			max += p.normal.y * this.min.y;

		}

		if(p.normal.z > 0) {

			min += p.normal.z * this.min.z;
			max += p.normal.z * this.max.z;

		} else {

			min += p.normal.z * this.max.z;
			max += p.normal.z * this.min.z;

		}

		return (min <= p.constant && max >= p.constant);

	};

	/**
		 * Checks if this box equals the given one.
		 *
		 * @param {Box3} v - A box.
		 * @return {Boolean} Whether this box equals the given one.
		 */

	Box3.prototype.equals = function equals (b) {

		return (b.min.equals(this.min) && b.max.equals(this.max));

	};

	/**
	 * A box.
	 *
	 * @type {Box3}
	 * @private
	 */

	var box = new Box3();

	/**
	 * A sphere.
	 */

	var Sphere = function Sphere(center, radius) {
		if ( center === void 0 ) center = new Vector3();
		if ( radius === void 0 ) radius = 0;


		/**
			 * The center.
			 *
			 * @type {Vector3}
			 */

		this.center = center;

		/**
			 * The radius.
			 *
			 * @type {Number}
			 */

		this.radius = radius;

	};

	/**
		 * Sets the center and the radius.
		 *
		 * @param {Vector3} center - The center.
		 * @param {Number} radius - The radius.
		 * @return {Sphere} This sphere.
		 */

	Sphere.prototype.set = function set (center, radius) {

		this.center.copy(center);
		this.radius = radius;

		return this;

	};

	/**
		 * Copies the given sphere.
		 *
		 * @param {Sphere} sphere - A sphere.
		 * @return {Sphere} This sphere.
		 */

	Sphere.prototype.copy = function copy (s) {

		this.center.copy(s.center);
		this.radius = s.radius;

		return this;

	};

	/**
		 * Clones this sphere.
		 *
		 * @return {Sphere} The cloned sphere.
		 */

	Sphere.prototype.clone = function clone () {

		return new this.constructor().copy(this);

	};

	/**
		 * Sets this sphere from points.
		 *
		 * @param {Vector3[]} points - The points.
		 * @param {Sphere} [target] - A target vector. If none is provided, a new one will be created.
		 * @return {Sphere} This sphere.
		 */

	Sphere.prototype.setFromPoints = function setFromPoints (points, target) {
			if ( target === void 0 ) target = box.setFromPoints(points).getCenter(this.center);


		var center = this.center;

		var maxRadiusSq = 0;
		var i, l;

		for(i = 0, l = points.length; i < l; ++i) {

			maxRadiusSq = Math.max(maxRadiusSq, center.distanceToSquared(points[i]));

		}

		this.radius = Math.sqrt(maxRadiusSq);

		return this;

	};

	/**
		 * Calculates the bounding box of this sphere.
		 *
		 * @param {Box3} [target] - A target sphere. If none is provided, a new one will be created.
		 * @return {Box3} The bounding box.
		 */

	Sphere.prototype.getBoundingBox = function getBoundingBox (target) {
			if ( target === void 0 ) target = new Box3();


		target.set(this.center, this.center);
		target.expandByScalar(this.radius);

		return target;

	};

	/**
		 * Checks if this sphere is empty.
		 *
		 * @return {Boolean} Whether this sphere is empty.
		 */

	Sphere.prototype.isEmpty = function isEmpty () {

		return (this.radius <= 0);

	};

	/**
		 * Translates this sphere.
		 *
		 * @param {Number} offset - An offset.
		 * @return {Sphere} This sphere.
		 */

	Sphere.prototype.translate = function translate (offset) {

		this.center.add(offset);

		return this;

	};

	/**
		 * Calculates the bounding box of this sphere.
		 *
		 * @param {Vector3} p - A point.
		 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
		 * @return {Vector3} The clamped point.
		 */

	Sphere.prototype.clampPoint = function clampPoint (p, target) {
			if ( target === void 0 ) target = new Vector3();


		var deltaLengthSq = this.center.distanceToSquared(p);

		target.copy(p);

		if(deltaLengthSq > (this.radius * this.radius)) {

			target.sub(this.center).normalize();
			target.multiplyScalar(this.radius).add(this.center);

		}

		return target;

	};

	/**
		 * Calculates the distance from this sphere to the given point.
		 *
		 * @param {Vector3} p - A point.
		 * @return {Number} The distance.
		 */

	Sphere.prototype.distanceToPoint = function distanceToPoint (p) {

		return (p.distanceTo(this.center) - this.radius);

	};

	/**
		 * Checks if the given point lies inside this sphere.
		 *
		 * @param {Vector3} p - A point.
		 * @return {Boolean} Whether this sphere contains the point.
		 */

	Sphere.prototype.containsPoint = function containsPoint (p) {

		return (p.distanceToSquared(this.center) <= (this.radius * this.radius));

	};

	/**
		 * Checks if the this sphere intersects with the given one.
		 *
		 * @param {Sphere} s - A sphere.
		 * @return {Boolean} Whether this sphere intersects with the given one.
		 */

	Sphere.prototype.intersectsSphere = function intersectsSphere (s) {

		var radiusSum = this.radius + s.radius;

		return s.center.distanceToSquared(this.center) <= (radiusSum * radiusSum);

	};

	/**
		 * Checks if the this sphere intersects with the given box.
		 *
		 * @param {Box3} b - A box.
		 * @return {Boolean} Whether this sphere intersects with the given box.
		 */

	Sphere.prototype.intersectsBox = function intersectsBox (b) {

		return b.intersectsSphere(this);

	};

	/**
		 * Checks if the this sphere intersects with the given plane.
		 *
		 * @param {Plane} p - A plane.
		 * @return {Boolean} Whether this sphere intersects with the given plane.
		 */

	Sphere.prototype.intersectsPlane = function intersectsPlane (p) {

		return (Math.abs(p.distanceToPoint(this.center)) <= this.radius);

	};

	/**
		 * Checks if this sphere equals the given one.
		 *
		 * @param {Sphere} s - A sphere.
		 * @return {Boolean} Whether the spheres are equal.
		 */

	Sphere.prototype.equals = function equals (s) {

		return (s.center.equals(this.center) && (s.radius === this.radius));

	};

	/**
	 * A vector with two components.
	 */

	var Vector2 = function Vector2(x, y) {
		if ( x === void 0 ) x = 0;
		if ( y === void 0 ) y = 0;


		/**
			 * The X component.
			 *
			 * @type {Number}
			 */

		this.x = x;

		/**
			 * The Y component.
			 *
			 * @type {Number}
			 */

		this.y = y;

	};

	var prototypeAccessors$1 = { width: { configurable: true },height: { configurable: true } };

	/**
		 * The width. This is an alias for X.
		 *
		 * @type {Number}
		 */

	prototypeAccessors$1.width.get = function () { return this.x; };

	/**
		 * Sets the width.
		 *
		 * @type {Number}
		 */

	prototypeAccessors$1.width.set = function (value) { return this.x = value; };

	/**
		 * The height. This is an alias for Y.
		 *
		 * @type {Number}
		 */

	prototypeAccessors$1.height.get = function () { return this.y; };

	/**
		 * Sets the height.
		 *
		 * @type {Number}
		 */

	prototypeAccessors$1.height.set = function (value) { return this.y = value; };

	/**
		 * Sets the values of this vector
		 *
		 * @param {Number} x - The X component.
		 * @param {Number} y - The Y component.
		 * @return {Vector2} This vector.
		 */

	Vector2.prototype.set = function set (x, y) {

		this.x = x;
		this.y = y;

		return this;

	};

	/**
		 * Copies the values of another vector.
		 *
		 * @param {Vector2} v - A vector.
		 * @return {Vector2} This vector.
		 */

	Vector2.prototype.copy = function copy (v) {

		this.x = v.x;
		this.y = v.y;

		return this;

	};

	/**
		 * Clones this vector.
		 *
		 * @return {Vector2} A clone of this vector.
		 */

	Vector2.prototype.clone = function clone () {

		return new this.constructor(this.x, this.y);

	};

	/**
		 * Copies values from an array.
		 *
		 * @param {Number[]} array - An array.
		 * @param {Number} offset - An offset.
		 * @return {Vector2} This vector.
		 */

	Vector2.prototype.fromArray = function fromArray (array, offset) {
			if ( offset === void 0 ) offset = 0;


		this.x = array[offset];
		this.y = array[offset + 1];

		return this;

	};

	/**
		 * Stores this vector in an array.
		 *
		 * @param {Array} [array] - A target array.
		 * @param {Number} offset - An offset.
		 * @return {Number[]} The array.
		 */

	Vector2.prototype.toArray = function toArray (array, offset) {
			if ( array === void 0 ) array = [];
			if ( offset === void 0 ) offset = 0;


		array[offset] = this.x;
		array[offset + 1] = this.y;

		return array;

	};

	/**
		 * Adds a vector to this one.
		 *
		 * @param {Vector2} v - The vector to add.
		 * @return {Vector2} This vector.
		 */

	Vector2.prototype.add = function add (v) {

		this.x += v.x;
		this.y += v.y;

		return this;

	};

	/**
		 * Adds a scalar to this vector.
		 *
		 * @param {Number} s - The scalar to add.
		 * @return {Vector2} This vector.
		 */

	Vector2.prototype.addScalar = function addScalar (s) {

		this.x += s;
		this.y += s;

		return this;

	};

	/**
		 * Sets this vector to the sum of two given vectors.
		 *
		 * @param {Vector2} a - A vector.
		 * @param {Vector2} b - Another vector.
		 * @return {Vector2} This vector.
		 */

	Vector2.prototype.addVectors = function addVectors (a, b) {

		this.x = a.x + b.x;
		this.y = a.y + b.y;

		return this;

	};

	/**
		 * Adds a scaled vector to this one.
		 *
		 * @param {Vector2} v - The vector to scale and add.
		 * @param {Number} s - A scalar.
		 * @return {Vector2} This vector.
		 */

	Vector2.prototype.addScaledVector = function addScaledVector (v, s) {

		this.x += v.x * s;
		this.y += v.y * s;

		return this;

	};

	/**
		 * Subtracts a vector from this vector.
		 *
		 * @param {Vector2} v - The vector to subtract.
		 * @return {Vector2} This vector.
		 */

	Vector2.prototype.sub = function sub (v) {

		this.x -= v.x;
		this.y -= v.y;

		return this;

	};

	/**
		 * Subtracts a scalar to this vector.
		 *
		 * @param {Number} s - The scalar to subtract.
		 * @return {Vector2} This vector.
		 */

	Vector2.prototype.subScalar = function subScalar (s) {

		this.x -= s;
		this.y -= s;

		return this;

	};

	/**
		 * Sets this vector to the difference between two given vectors.
		 *
		 * @param {Vector2} a - A vector.
		 * @param {Vector2} b - A second vector.
		 * @return {Vector2} This vector.
		 */

	Vector2.prototype.subVectors = function subVectors (a, b) {

		this.x = a.x - b.x;
		this.y = a.y - b.y;

		return this;

	};

	/**
		 * Multiplies this vector with another vector.
		 *
		 * @param {Vector2} v - A vector.
		 * @return {Vector2} This vector.
		 */

	Vector2.prototype.multiply = function multiply (v) {

		this.x *= v.x;
		this.y *= v.y;

		return this;

	};

	/**
		 * Multiplies this vector with a given scalar.
		 *
		 * @param {Number} s - A scalar.
		 * @return {Vector2} This vector.
		 */

	Vector2.prototype.multiplyScalar = function multiplyScalar (s) {

		this.x *= s;
		this.y *= s;

		return this;

	};

	/**
		 * Divides this vector by another vector.
		 *
		 * @param {Vector2} v - A vector.
		 * @return {Vector2} This vector.
		 */

	Vector2.prototype.divide = function divide (v) {

		this.x /= v.x;
		this.y /= v.y;

		return this;

	};

	/**
		 * Divides this vector by a given scalar.
		 *
		 * @param {Number} s - A scalar.
		 * @return {Vector2} This vector.
		 */

	Vector2.prototype.divideScalar = function divideScalar (s) {

		this.x /= s;
		this.y /= s;

		return this;

	};

	/**
		 * Applies the given matrix to this vector.
		 *
		 * @param {Matrix3} m - A matrix.
		 * @return {Vector2} This vector.
		 */

	Vector2.prototype.applyMatrix3 = function applyMatrix3 (m) {

		var x = this.x, y = this.y;
		var e = m.elements;

		this.x = e[0] * x + e[3] * y + e[6];
		this.y = e[1] * x + e[4] * y + e[7];

		return this;

	};

	/**
		 * Calculates the dot product with another vector.
		 *
		 * @param {Vector2} v - A vector.
		 * @return {Number} The dot product.
		 */

	Vector2.prototype.dot = function dot (v) {

		return this.x * v.x + this.y * v.y;

	};

	/**
		 * Calculates the Manhattan length of this vector.
		 *
		 * @return {Number} The length.
		 */

	Vector2.prototype.lengthManhattan = function lengthManhattan () {

		return Math.abs(this.x) + Math.abs(this.y);

	};

	/**
		 * Calculates the squared length of this vector.
		 *
		 * @return {Number} The squared length.
		 */

	Vector2.prototype.lengthSquared = function lengthSquared () {

		return this.x * this.x + this.y * this.y;

	};

	/**
		 * Calculates the length of this vector.
		 *
		 * @return {Number} The length.
		 */

	Vector2.prototype.length = function length () {

		return Math.sqrt(this.x * this.x + this.y * this.y);

	};

	/**
		 * Calculates the Manhattan distance to a given vector.
		 *
		 * @param {Vector2} v - A vector.
		 * @return {Number} The squared distance.
		 */

	Vector2.prototype.distanceToManhattan = function distanceToManhattan (v) {

		return Math.abs(this.x - v.x) + Math.abs(this.y - v.y);

	};

	/**
		 * Calculates the squared distance to a given vector.
		 *
		 * @param {Vector2} v - A vector.
		 * @return {Number} The squared distance.
		 */

	Vector2.prototype.distanceToSquared = function distanceToSquared (v) {

		var dx = this.x - v.x;
		var dy = this.y - v.y;

		return dx * dx + dy * dy;

	};

	/**
		 * Calculates the distance to a given vector.
		 *
		 * @param {Vector2} v - A vector.
		 * @return {Number} The distance.
		 */

	Vector2.prototype.distanceTo = function distanceTo (v) {

		return Math.sqrt(this.distanceToSquared(v));

	};

	/**
		 * Normalizes this vector.
		 *
		 * @return {Vector2} This vector.
		 */

	Vector2.prototype.normalize = function normalize () {

		return this.divideScalar(this.length());

	};

	/**
		 * Sets the length of this vector.
		 *
		 * @param {Number} length - The new length.
		 * @return {Vector2} This vector.
		 */

	Vector2.prototype.setLength = function setLength (length) {

		return this.normalize().multiplyScalar(length);

	};

	/**
		 * Adopts the min value for each component of this vector and the given one.
		 *
		 * @param {Vector2} v - A vector.
		 * @return {Vector2} This vector.
		 */

	Vector2.prototype.min = function min (v) {

		this.x = Math.min(this.x, v.x);
		this.y = Math.min(this.y, v.y);

		return this;

	};

	/**
		 * adopts the max value for each component of this vector and the given one.
		 *
		 * @param {Vector2} v - A vector.
		 * @return {Vector2} This vector.
		 */

	Vector2.prototype.max = function max (v) {

		this.x = Math.max(this.x, v.x);
		this.y = Math.max(this.y, v.y);

		return this;

	};

	/**
		 * Clamps this vector.
		 *
		 * @param {Vector2} min - A vector, assumed to be smaller than max.
		 * @param {Vector2} max - A vector, assumed to be greater than min.
		 * @return {Vector2} This vector.
		 */

	Vector2.prototype.clamp = function clamp (min, max) {

		this.x = Math.max(min.x, Math.min(max.x, this.x));
		this.y = Math.max(min.y, Math.min(max.y, this.y));

		return this;

	};

	/**
		 * Floors this vector.
		 *
		 * @return {Vector2} This vector.
		 */

	Vector2.prototype.floor = function floor () {

		this.x = Math.floor(this.x);
		this.y = Math.floor(this.y);

		return this;

	};

	/**
		 * Ceils this vector.
		 *
		 * @return {Vector2} This vector.
		 */

	Vector2.prototype.ceil = function ceil () {

		this.x = Math.ceil(this.x);
		this.y = Math.ceil(this.y);

		return this;

	};

	/**
		 * Rounds this vector.
		 *
		 * @return {Vector2} This vector.
		 */

	Vector2.prototype.round = function round () {

		this.x = Math.round(this.x);
		this.y = Math.round(this.y);

		return this;

	};

	/**
		 * Negates this vector.
		 *
		 * @return {Vector2} This vector.
		 */

	Vector2.prototype.negate = function negate () {

		this.x = -this.x;
		this.y = -this.y;

		return this;

	};

	/**
		 * Computes the angle in radians with respect to the positive X-axis.
		 *
		 * @return {Number} The angle.
		 */

	Vector2.prototype.angle = function angle () {

		var angle = Math.atan2(this.y, this.x);

		if(angle < 0) { angle += 2 * Math.PI; }

		return angle;

	};

	/**
		 * Lerps towards the given vector.
		 *
		 * @param {Vector2} v - The target vector.
		 * @param {Number} alpha - The lerp factor.
		 * @return {Vector2} This vector.
		 */

	Vector2.prototype.lerp = function lerp (v, alpha) {

		this.x += (v.x - this.x) * alpha;
		this.y += (v.y - this.y) * alpha;

		return this;

	};

	/**
		 * Sets this vector to the lerp result of the given vectors.
		 *
		 * @param {Vector2} v1 - A base vector.
		 * @param {Vector2} v2 - The target vector.
		 * @param {Number} alpha - The lerp factor.
		 * @return {Vector2} This vector.
		 */

	Vector2.prototype.lerpVectors = function lerpVectors (v1, v2, alpha) {

		return this.subVectors(v2, v1).multiplyScalar(alpha).add(v1);

	};

	/**
		 * Rotates this vector around a given center.
		 *
		 * @param {Vector2} center - The center.
		 * @param {Number} angle - The rotation in radians.
		 * @return {Vector2} This vector.
		 */

	Vector2.prototype.rotateAround = function rotateAround (center, angle) {

		var c = Math.cos(angle), s = Math.sin(angle);

		var x = this.x - center.x;
		var y = this.y - center.y;

		this.x = x * c - y * s + center.x;
		this.y = x * s + y * c + center.y;

		return this;

	};

	/**
		 * Checks if this vector equals the given one.
		 *
		 * @param {Vector2} v - A vector.
		 * @return {Boolean} Whether this vector equals the given one.
		 */

	Vector2.prototype.equals = function equals (v) {

		return (v.x === this.x && v.y === this.y);

	};

	Object.defineProperties( Vector2.prototype, prototypeAccessors$1 );

	/**
	 * A vector.
	 *
	 * @type {Vector2}
	 * @private
	 */

	var v = new Vector2();

	/**
	 * A 2D box.
	 */

	var Box2 = function Box2(
		min,
		max
	) {
		if ( min === void 0 ) min = new Vector2(Infinity, Infinity);
		if ( max === void 0 ) max = new Vector2(-Infinity, -Infinity);


		/**
			 * The lower bounds.
			 *
			 * @type {Vector2}
			 */

		this.min = min;

		/**
			 * The upper bounds.
			 *
			 * @type {Vector2}
			 */

		this.max = max;

	};

	/**
		 * Sets the values of this box.
		 *
		 * @param {Vector2} min - The lower bounds.
		 * @param {Vector2} max - The upper bounds.
		 * @return {Box2} This box.
		 */

	Box2.prototype.set = function set (min, max) {

		this.min.copy(min);
		this.max.copy(max);

		return this;

	};

	/**
		 * Copies the values of a given box.
		 *
		 * @param {Box2} b - A box.
		 * @return {Box2} This box.
		 */

	Box2.prototype.copy = function copy (b) {

		this.min.copy(b.min);
		this.max.copy(b.max);

		return this;

	};

	/**
		 * Clones this box.
		 *
		 * @return {Box2} A clone of this box.
		 */

	Box2.prototype.clone = function clone () {

		return new this.constructor().copy(this);

	};

	/**
		 * Makes this box empty.
		 *
		 * The lower bounds are set to infinity and the upper bounds to negative
		 * infinity to create an infinitely small box.
		 *
		 * @return {Box2} This box.
		 */

	Box2.prototype.makeEmpty = function makeEmpty () {

		this.min.x = this.min.y = Infinity;
		this.max.x = this.max.y = -Infinity;

		return this;

	};

	/**
		 * Indicates whether this box is truly empty.
		 *
		 * This is a more robust check for emptiness since the volume can get positive
		 * with two negative axes.
		 *
		 * @return {Box2} This box.
		 */

	Box2.prototype.isEmpty = function isEmpty () {

		return (
			this.max.x < this.min.x ||
			this.max.y < this.min.y
		);

	};

	/**
		 * Computes the center of this box.
		 *
		 * @param {Vector2} [target] - A target vector. If none is provided, a new one will be created.
		 * @return {Vector2} A vector that describes the center of this box.
		 */

	Box2.prototype.getCenter = function getCenter (target) {
			if ( target === void 0 ) target = new Vector2();


		return !this.isEmpty() ?
			target.addVectors(this.min, this.max).multiplyScalar(0.5) :
			target.set(0, 0);

	};

	/**
		 * Computes the size of this box.
		 *
		 * @param {Vector2} [target] - A target vector. If none is provided, a new one will be created.
		 * @return {Vector2} A vector that describes the size of this box.
		 */

	Box2.prototype.getSize = function getSize (target) {
			if ( target === void 0 ) target = new Vector2();


		return !this.isEmpty() ?
			target.subVectors(this.max, this.min) :
			target.set(0, 0);

	};

	/**
		 * Computes the bounding sphere of this box.
		 *
		 * @param {Sphere} [target] - A target sphere. If none is provided, a new one will be created.
		 * @return {Sphere} The bounding sphere of this box.
		 */

	Box2.prototype.getBoundingSphere = function getBoundingSphere (target) {
			if ( target === void 0 ) target = new Sphere();


		this.getCenter(target.center);

		target.radius = this.getSize(v).length() * 0.5;

		return target;

	};

	/**
		 * Expands this box by the given point.
		 *
		 * @param {Vector2} p - A point.
		 * @return {Box2} This box.
		 */

	Box2.prototype.expandByPoint = function expandByPoint (p) {

		this.min.min(p);
		this.max.max(p);

		return this;

	};

	/**
		 * Expands this box by the given vector.
		 *
		 * @param {Vector2} v - A vector.
		 * @return {Box2} This box.
		 */

	Box2.prototype.expandByVector = function expandByVector (v) {

		this.min.sub(v);
		this.max.add(v);

		return this;

	};

	/**
		 * Expands this box by the given scalar.
		 *
		 * @param {Number} s - A scalar.
		 * @return {Box2} This box.
		 */

	Box2.prototype.expandByScalar = function expandByScalar (s) {

		this.min.addScalar(-s);
		this.max.addScalar(s);

		return this;

	};

	/**
		 * Defines this box by the given points.
		 *
		 * @param {Vector2[]} points - The points.
		 * @return {Box2} This box.
		 */

	Box2.prototype.setFromPoints = function setFromPoints (points) {
			var this$1 = this;


		var i, l;

		this.min.set(0, 0);
		this.max.set(0, 0);

		for(i = 0, l = points.length; i < l; ++i) {

			this$1.expandByPoint(points[i]);

		}

		return this;

	};

	/**
		 * Defines this box by the given center and size.
		 *
		 * @param {Vector2} center - The center.
		 * @param {Number} size - The size.
		 * @return {Box2} This box.
		 */

	Box2.prototype.setFromCenterAndSize = function setFromCenterAndSize (center, size) {

		var halfSize = v.copy(size).multiplyScalar(0.5);

		this.min.copy(center).sub(halfSize);
		this.max.copy(center).add(halfSize);

		return this;

	};

	/**
		 * Clamps the given point to the boundaries of this box.
		 *
		 * @param {Vector2} p - A point.
		 * @param {Vector2} [target] - A target vector. If none is provided, a new one will be created.
		 * @return {Vector2} The clamped point.
		 */

	Box2.prototype.clampPoint = function clampPoint (point, target) {
			if ( target === void 0 ) target = new Vector2();


		return target.copy(point).clamp(this.min, this.max);

	};

	/**
		 * Calculates the distance from this box to the given point.
		 *
		 * @param {Vector2} p - A point.
		 * @return {Number} The distance.
		 */

	Box2.prototype.distanceToPoint = function distanceToPoint (p) {

		var clampedPoint = v.copy(p).clamp(this.min, this.max);

		return clampedPoint.sub(p).length();

	};

	/**
		 * Translates this box.
		 *
		 * @param {Vector2} offset - The offset.
		 * @return {Box2} This box.
		 */

	Box2.prototype.translate = function translate (offset) {

		this.min.add(offset);
		this.max.add(offset);

		return this;

	};

	/**
		 * Expands this box by combining it with the given one.
		 *
		 * @param {Box2} b - A box.
		 * @return {Box2} This box.
		 */

	Box2.prototype.intersect = function intersect (b) {

		this.min.max(b.min);
		this.max.min(b.max);

		/* Ensure that if there is no overlap, the result is fully empty to prevent
		subsequent intersections to erroneously return valid values. */
		if(this.isEmpty()) { this.makeEmpty(); }

		return this;

	};

	/**
		 * Expands this box by combining it with the given one.
		 *
		 * @param {Box2} b - A box.
		 * @return {Box2} This box.
		 */

	Box2.prototype.union = function union (b) {

		this.min.min(b.min);
		this.max.max(b.max);

		return this;

	};

	/**
		 * Checks if the given point lies inside this box.
		 *
		 * @param {Vector2} p - A point.
		 * @return {Boolean} Whether this box contains the point.
		 */

	Box2.prototype.containsPoint = function containsPoint (p) {

		return !(
			p.x < this.min.x || p.x > this.max.x ||
			p.y < this.min.y || p.y > this.max.y
		);

	};

	/**
		 * Checks if the given box lies inside this box.
		 *
		 * @param {Vector2} b - A box.
		 * @return {Boolean} Whether this box contains the given one.
		 */

	Box2.prototype.containsBox = function containsBox (b) {

		return (
			this.min.x <= b.min.x && b.max.x <= this.max.x &&
			this.min.y <= b.min.y && b.max.y <= this.max.y
		);

	};

	/**
		 * Checks if this box intersects with the given one.
		 *
		 * @param {Box2} b - A box.
		 * @return {Boolean} Whether the boxes intersect.
		 */

	Box2.prototype.intersectsBox = function intersectsBox (b) {

		return !(
			b.max.x < this.min.x || b.min.x > this.max.x ||
			b.max.y < this.min.y || b.min.y > this.max.y
		);

	};

	/**
		 * Checks if this box equals the given one.
		 *
		 * @param {Box2} v - A box.
		 * @return {Boolean} Whether this box equals the given one.
		 */

	Box2.prototype.equals = function equals (b) {

		return (b.min.equals(this.min) && b.max.equals(this.max));

	};

	/**
	 * A cylindrical coordinate system.
	 *
	 * For details see: https://en.wikipedia.org/wiki/Cylindrical_coordinate_system
	 */

	var Cylindrical = function Cylindrical(radius, theta, y) {
		if ( radius === void 0 ) radius = 1;
		if ( theta === void 0 ) theta = 0;
		if ( y === void 0 ) y = 0;


		/**
			 * The distance from the origin to a point in the XZ-plane.
			 *
			 * @type {Number}
			 */

		this.radius = radius;

		/**
			 * The counterclockwise angle in the XZ-plane measured in radians from the
			 * positive Z-axis.
			 *
			 * @type {Number}
			 */

		this.theta = theta;

		/**
			 * The height above the XZ-plane.
			 *
			 * @type {Number}
			 */

		this.y = y;

	};

	/**
		 * Sets the values of this cylindrical system.
		 *
		 * @param {Number} radius - The radius.
		 * @param {Number} theta - Theta.
		 * @param {Number} y - The height.
		 * @return {Cylindrical} This cylindrical system.
		 */

	Cylindrical.prototype.set = function set (radius, theta, y) {

		this.radius = radius;
		this.theta = theta;
		this.y = y;

		return this;

	};

	/**
		 * Copies the values of the given cylindrical system.
		 *
		 * @param {Cylindrical} c - A cylindrical system.
		 * @return {Cylindrical} This cylindrical system.
		 */

	Cylindrical.prototype.copy = function copy (c) {

		this.radius = c.radius;
		this.theta = c.theta;
		this.y = c.y;

		return this;

	};

	/**
		 * Clones this cylindrical system.
		 *
		 * @return {Cylindrical} The cloned cylindrical system.
		 */

	Cylindrical.prototype.clone = function clone () {

		return new this.constructor().copy(this);

	};

	/**
		 * Sets the values of this cylindrical system.
		 *
		 * @param {Vector3} v - The vector
		 * @return {Cylindrical} This cylindrical system.
		 */

	Cylindrical.prototype.setFromVector3 = function setFromVector3 (v) {

		this.radius = Math.sqrt(v.x * v.x + v.z * v.z);
		this.theta = Math.atan2(v.x, v.z);
		this.y = v.y;

		return this;

	};

	/**
	 * A 3x3 matrix.
	 */

	var Matrix3 = function Matrix3() {

		/**
			 * The matrix elements.
			 *
			 * @type {Float32Array}
			 */

		this.elements = new Float32Array([

			1, 0, 0,
			0, 1, 0,
			0, 0, 1

		]);

	};

	/**
		 * Sets the values of this matrix.
		 *
		 * @param {Number} m00 - The value of the first row, first column.
		 * @param {Number} m01 - The value of the first row, second column.
		 * @param {Number} m02 - The value of the first row, third column.
		 * @param {Number} m10 - The value of the second row, first column.
		 * @param {Number} m11 - The value of the second row, second column.
		 * @param {Number} m12 - The value of the second row, third column.
		 * @param {Number} m20 - The value of the third row, first column.
		 * @param {Number} m21 - The value of the third row, second column.
		 * @param {Number} m22 - The value of the third row, third column.
		 * @return {Matrix3} This matrix.
		 */

	Matrix3.prototype.set = function set (m00, m01, m02, m10, m11, m12, m20, m21, m22) {

		var te = this.elements;

		te[0] = m00; te[3] = m01; te[6] = m02;
		te[1] = m10; te[4] = m11; te[7] = m12;
		te[2] = m20; te[5] = m21; te[8] = m22;

		return this;

	};

	/**
		 * Sets this matrix to the identity matrix.
		 *
		 * @return {Matrix3} This matrix.
		 */

	Matrix3.prototype.identity = function identity () {

		this.set(

			1, 0, 0,
			0, 1, 0,
			0, 0, 1

		);

		return this;

	};

	/**
		 * Copies the values of a given matrix.
		 *
		 * @param {Matrix3} matrix - A matrix.
		 * @return {Matrix3} This matrix.
		 */

	Matrix3.prototype.copy = function copy (matrix) {

		var me = matrix.elements;
		var te = this.elements;

		te[0] = me[0]; te[1] = me[1]; te[2] = me[2];
		te[3] = me[3]; te[4] = me[4]; te[5] = me[5];
		te[6] = me[6]; te[7] = me[7]; te[8] = me[8];

		return this;

	};

	/**
		 * Clones this matrix.
		 *
		 * @return {Matrix3} A clone of this matrix.
		 */

	Matrix3.prototype.clone = function clone () {

		return new this.constructor().fromArray(this.elements);

	};

	/**
		 * Copies the values of a given array.
		 *
		 * @param {Number[]} array - An array.
		 * @param {Number} [offset=0] - An offset into the array.
		 * @return {Matrix3} This matrix.
		 */

	Matrix3.prototype.fromArray = function fromArray (array, offset) {
			if ( offset === void 0 ) offset = 0;


		var te = this.elements;

		var i;

		for(i = 0; i < 9; ++i) {

			te[i] = array[i + offset];

		}

		return this;

	};

	/**
		 * Stores this matrix in an array.
		 *
		 * @param {Number[]} [array] - A target array.
		 * @param {Number} [offset=0] - An offset into the array.
		 * @return {Number[]} The array.
		 */

	Matrix3.prototype.toArray = function toArray (array, offset) {
			if ( array === void 0 ) array = [];
			if ( offset === void 0 ) offset = 0;


		var te = this.elements;

		var i;

		for(i = 0; i < 9; ++i) {

			array[i + offset] = te[i];

		}

		return array;

	};

	/**
		 * Sets this matrix to the product of the given matrices.
		 *
		 * @param {Matrix3} a - A matrix.
		 * @param {Matrix3} b - A matrix.
		 * @return {Matrix3} This matrix.
		 */

	Matrix3.prototype.multiplyMatrices = function multiplyMatrices (a, b) {

		var ae = a.elements;
		var be = b.elements;
		var te = this.elements;

		var a11 = ae[0], a12 = ae[3], a13 = ae[6];
		var a21 = ae[1], a22 = ae[4], a23 = ae[7];
		var a31 = ae[2], a32 = ae[5], a33 = ae[8];

		var b11 = be[0], b12 = be[3], b13 = be[6];
		var b21 = be[1], b22 = be[4], b23 = be[7];
		var b31 = be[2], b32 = be[5], b33 = be[8];

		te[0] = a11 * b11 + a12 * b21 + a13 * b31;
		te[3] = a11 * b12 + a12 * b22 + a13 * b32;
		te[6] = a11 * b13 + a12 * b23 + a13 * b33;

		te[1] = a21 * b11 + a22 * b21 + a23 * b31;
		te[4] = a21 * b12 + a22 * b22 + a23 * b32;
		te[7] = a21 * b13 + a22 * b23 + a23 * b33;

		te[2] = a31 * b11 + a32 * b21 + a33 * b31;
		te[5] = a31 * b12 + a32 * b22 + a33 * b32;
		te[8] = a31 * b13 + a32 * b23 + a33 * b33;

		return this;

	};

	/**
		 * Multiplies this matrix with a given one.
		 *
		 * @param {Matrix3} m - A matrix.
		 * @return {Matrix3} This matrix.
		 */

	Matrix3.prototype.multiply = function multiply (m) {

		return this.multiplyMatrices(this, m);

	};

	/**
		 * Multiplies a given matrix with this one.
		 *
		 * @param {Matrix3} m - A matrix.
		 * @return {Matrix3} This matrix.
		 */

	Matrix3.prototype.premultiply = function premultiply (m) {

		return this.multiplyMatrices(m, this);

	};

	/**
		 * Multiplies this matrix with a given scalar.
		 *
		 * @param {Number} m - A scalar.
		 * @return {Matrix3} This matrix.
		 */

	Matrix3.prototype.multiplyScalar = function multiplyScalar (s) {

		var te = this.elements;

		te[0] *= s; te[3] *= s; te[6] *= s;
		te[1] *= s; te[4] *= s; te[7] *= s;
		te[2] *= s; te[5] *= s; te[8] *= s;

		return this;

	};

	/**
		 * Calculates the determinant of this matrix.
		 *
		 * @return {Number} The determinant.
		 */

	Matrix3.prototype.determinant = function determinant () {

		var te = this.elements;

		var a = te[0], b = te[1], c = te[2];
		var d = te[3], e = te[4], f = te[5];
		var g = te[6], h = te[7], i = te[8];

		return (

			a * e * i -
			a * f * h -
			b * d * i +
			b * f * g +
			c * d * h -
			c * e * g

		);

	};

	/**
		 * Inverts the given matrix and stores the result in this matrix.
		 *
		 * @param {Matrix3} matrix - The matrix that should be inverted.
		 * @return {Matrix3} This matrix.
		 */

	Matrix3.prototype.getInverse = function getInverse (matrix) {

		var me = matrix.elements;
		var te = this.elements;

		var n11 = me[0], n21 = me[1], n31 = me[2];
		var n12 = me[3], n22 = me[4], n32 = me[5];
		var n13 = me[6], n23 = me[7], n33 = me[8];

		var t11 = n33 * n22 - n32 * n23;
		var t12 = n32 * n13 - n33 * n12;
		var t13 = n23 * n12 - n22 * n13;

		var det = n11 * t11 + n21 * t12 + n31 * t13;

		var invDet;

		if(det !== 0) {

			invDet = 1.0 / det;

			te[0] = t11 * invDet;
			te[1] = (n31 * n23 - n33 * n21) * invDet;
			te[2] = (n32 * n21 - n31 * n22) * invDet;

			te[3] = t12 * invDet;
			te[4] = (n33 * n11 - n31 * n13) * invDet;
			te[5] = (n31 * n12 - n32 * n11) * invDet;

			te[6] = t13 * invDet;
			te[7] = (n21 * n13 - n23 * n11) * invDet;
			te[8] = (n22 * n11 - n21 * n12) * invDet;

		} else {

			console.error("Can't invert matrix, determinant is zero", matrix);

			this.identity();

		}

		return this;

	};

	/**
		 * Transposes this matrix.
		 *
		 * @return {Matrix3} This matrix.
		 */

	Matrix3.prototype.transpose = function transpose () {

		var me = this.elements;

		var t;

		t = me[1]; me[1] = me[3]; me[3] = t;
		t = me[2]; me[2] = me[6]; me[6] = t;
		t = me[5]; me[5] = me[7]; me[7] = t;

		return this;

	};

	/**
		 * Scales this matrix.
		 *
		 * @param {Number} sx - The X scale.
		 * @param {Number} sy - The Y scale.
		 * @return {Matrix3} This matrix.
		 */

	Matrix3.prototype.scale = function scale (sx, sy) {

		var te = this.elements;

		te[0] *= sx; te[3] *= sx; te[6] *= sx;
		te[1] *= sy; te[4] *= sy; te[7] *= sy;

		return this;

	};

	/**
		 * Rotates this matrix.
		 *
		 * @param {Number} theta - The rotation.
		 * @return {Matrix3} This matrix.
		 */

	Matrix3.prototype.rotate = function rotate (theta) {

		var c = Math.cos(theta);
		var s = Math.sin(theta);

		var te = this.elements;

		var a11 = te[0], a12 = te[3], a13 = te[6];
		var a21 = te[1], a22 = te[4], a23 = te[7];

		te[0] = c * a11 + s * a21;
		te[3] = c * a12 + s * a22;
		te[6] = c * a13 + s * a23;

		te[1] = -s * a11 + c * a21;
		te[4] = -s * a12 + c * a22;
		te[7] = -s * a13 + c * a23;

		return this;

	};

	/**
		 * Translates this matrix.
		 *
		 * @param {Number} tx - The X offset.
		 * @param {Number} ty - The Y offset.
		 * @return {Matrix3} This matrix.
		 */

	Matrix3.prototype.translate = function translate (tx, ty) {

		var te = this.elements;

		te[0] += tx * te[2]; te[3] += tx * te[5]; te[6] += tx * te[8];
		te[1] += ty * te[2]; te[4] += ty * te[5]; te[7] += ty * te[8];

		return this;

	};

	/**
		 * Checks if this matrix equals the given one.
		 *
		 * @param {Matrix3} m - A matrix.
		 * @return {Boolean} Whether the matrix are equal.
		 */

	Matrix3.prototype.equals = function equals (matrix) {

		var te = this.elements;
		var me = matrix.elements;

		var result = true;
		var i;

		for(i = 0; result && i < 9; ++i) {

			if(te[i] !== me[i]) {

				result = false;

			}

		}

		return result;

	};

	/**
	 * An enumeration of Euler rotation orders.
	 *
	 * @type {Object}
	 * @property {String} XYZ - X -> Y -> Z.
	 * @property {String} YZX - Y -> Z -> X.
	 * @property {String} ZXY - Z -> X -> Y.
	 * @property {String} XZY - X -> Z -> Y.
	 * @property {String} YXZ - Y -> X -> Z.
	 * @property {String} ZYX - Z -> Y -> X.
	 */

	var RotationOrder = {

		XYZ: "XYZ",
		YZX: "YZX",
		ZXY: "ZXY",
		XZY: "XZY",
		YXZ: "YXZ",
		ZYX: "ZYX"

	};

	/**
	 * A vector.
	 *
	 * @type {Vector3}
	 * @private
	 */

	var v$2 = new Vector3();

	/**
	 * A quaternion.
	 */

	var Quaternion = function Quaternion(x, y, z, w) {
		if ( x === void 0 ) x = 0;
		if ( y === void 0 ) y = 0;
		if ( z === void 0 ) z = 0;
		if ( w === void 0 ) w = 0;


		/**
			 * The X component.
			 *
			 * @type {Number}
			 */

		this.x = x;

		/**
			 * The X component.
			 *
			 * @type {Number}
			 */

		this.y = x;

		/**
			 * The X component.
			 *
			 * @type {Number}
			 */

		this.z = x;

		/**
			 * The X component.
			 *
			 * @type {Number}
			 */

		this.w = x;

	};

	/**
		 * Sets the components of this quaternion.
		 *
		 * @param {Number} x - The X component.
		 * @param {Number} y - The X component.
		 * @param {Number} z - The X component.
		 * @param {Number} w - The X component.
		 * @return {Quaternion} This quaternion.
		 */

	Quaternion.prototype.set = function set (x, y, z, w) {

		this.x = x;
		this.y = y;
		this.z = z;
		this.w = w;

		return this;

	};

	/**
		 * Copies the components of the given quaternion.
		 *
		 * @param {Quaternion} q - The quaternion.
		 * @return {Quaternion} This quaternion.
		 */

	Quaternion.prototype.copy = function copy (q) {

		this.x = q.x;
		this.y = q.y;
		this.z = q.z;
		this.w = q.w;

		return this;

	};

	/**
		 * Clones this quaternion.
		 *
		 * @return {Quaternion} The cloned quaternion.
		 */

	Quaternion.prototype.clone = function clone () {

		return new this.constructor(this.x, this.y, this.z, this.w);

	};

	/**
		 * Copies values from an array.
		 *
		 * @param {Number[]} array - An array.
		 * @param {Number} offset - An offset.
		 * @return {Quaternion} This quaternion.
		 */

	Quaternion.prototype.fromArray = function fromArray (array, offset) {
			if ( offset === void 0 ) offset = 0;


		this.x = array[offset];
		this.y = array[offset + 1];
		this.z = array[offset + 2];
		this.w = array[offset + 3];

		return this;

	};

	/**
		 * Stores this quaternion in an array.
		 *
		 * @param {Array} [array] - A target array.
		 * @param {Number} offset - An offset.
		 * @return {Number[]} The array.
		 */

	Quaternion.prototype.toArray = function toArray (array, offset) {
			if ( array === void 0 ) array = [];
			if ( offset === void 0 ) offset = 0;


		array[offset] = this.x;
		array[offset + 1] = this.y;
		array[offset + 2] = this.z;
		array[offset + 3] = this.w;

		return array;

	};

	/**
		 * Sets the components of this quaternion based on the given Euler angles.
		 *
		 * For more details see: https://goo.gl/XRD1kr
		 *
		 * @return {Euler} euler - The euler angles.
		 * @return {Quaternion} This quaternion.
		 */

	Quaternion.prototype.setFromEuler = function setFromEuler (euler) {

		var x = euler.x;
		var y = euler.y;
		var z = euler.z;

		var cos = Math.cos;
		var sin = Math.sin;

		var c1 = cos(x / 2);
		var c2 = cos(y / 2);
		var c3 = cos(z / 2);

		var s1 = sin(x / 2);
		var s2 = sin(y / 2);
		var s3 = sin(z / 2);

		switch(euler.order) {

			case RotationOrder.XYZ:
				this.x = s1 * c2 * c3 + c1 * s2 * s3;
				this.y = c1 * s2 * c3 - s1 * c2 * s3;
				this.z = c1 * c2 * s3 + s1 * s2 * c3;
				this.w = c1 * c2 * c3 - s1 * s2 * s3;
				break;

			case RotationOrder.YXZ:
				this.x = s1 * c2 * c3 + c1 * s2 * s3;
				this.y = c1 * s2 * c3 - s1 * c2 * s3;
				this.z = c1 * c2 * s3 - s1 * s2 * c3;
				this.w = c1 * c2 * c3 + s1 * s2 * s3;
				break;

			case RotationOrder.ZXY:
				this.x = s1 * c2 * c3 - c1 * s2 * s3;
				this.y = c1 * s2 * c3 + s1 * c2 * s3;
				this.z = c1 * c2 * s3 + s1 * s2 * c3;
				this.w = c1 * c2 * c3 - s1 * s2 * s3;
				break;

			case RotationOrder.ZYX:
				this.x = s1 * c2 * c3 - c1 * s2 * s3;
				this.y = c1 * s2 * c3 + s1 * c2 * s3;
				this.z = c1 * c2 * s3 - s1 * s2 * c3;
				this.w = c1 * c2 * c3 + s1 * s2 * s3;
				break;

			case RotationOrder.YZX:
				this.x = s1 * c2 * c3 + c1 * s2 * s3;
				this.y = c1 * s2 * c3 + s1 * c2 * s3;
				this.z = c1 * c2 * s3 - s1 * s2 * c3;
				this.w = c1 * c2 * c3 - s1 * s2 * s3;
				break;

			case RotationOrder.XZY:
				this.x = s1 * c2 * c3 - c1 * s2 * s3;
				this.y = c1 * s2 * c3 - s1 * c2 * s3;
				this.z = c1 * c2 * s3 + s1 * s2 * c3;
				this.w = c1 * c2 * c3 + s1 * s2 * s3;
				break;

		}

		return this;

	};

	/**
		 * Sets the components of this quaternion based on a given axis angle.
		 *
		 * For more information see:
		 *  http://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToQuaternion/index.htm
		 *
		 * @param {Vector3} axis - The axis. Assumed to be normalized.
		 * @param {Number} angle - The angle in radians.
		 * @return {Quaternion} This quaternion.
		 */

	Quaternion.prototype.setFromAxisAngle = function setFromAxisAngle (axis, angle) {

		var halfAngle = angle / 2.0;
		var s = Math.sin(halfAngle);

		this.x = axis.x * s;
		this.y = axis.y * s;
		this.z = axis.z * s;
		this.w = Math.cos(halfAngle);

		return this;

	};

	/**
		 * Sets the components of this quaternion based on a given rotation matrix.
		 *
		 * For more information see:
		 *  http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm
		 *
		 * @param {Matrix4} m - The rotation matrix. The upper 3x3 is assumed to be a pure rotation matrix (i.e. unscaled).
		 * @return {Quaternion} This quaternion.
		 */

	Quaternion.prototype.setFromRotationMatrix = function setFromRotationMatrix (m) {

		var te = m.elements;

		var m00 = te[0], m01 = te[4], m02 = te[8];
		var m10 = te[1], m11 = te[5], m12 = te[9];
		var m20 = te[2], m21 = te[6], m22 = te[10];

		var trace = m00 + m11 + m22;

		var s;

		if(trace > 0) {

			s = 0.5 / Math.sqrt(trace + 1.0);

			this.w = 0.25 / s;
			this.x = (m21 - m12) * s;
			this.y = (m02 - m20) * s;
			this.z = (m10 - m01) * s;

		} else if(m00 > m11 && m00 > m22) {

			s = 2.0 * Math.sqrt(1.0 + m00 - m11 - m22);

			this.w = (m21 - m12) / s;
			this.x = 0.25 * s;
			this.y = (m01 + m10) / s;
			this.z = (m02 + m20) / s;

		} else if(m11 > m22) {

			s = 2.0 * Math.sqrt(1.0 + m11 - m00 - m22);

			this.w = (m02 - m20) / s;
			this.x = (m01 + m10) / s;
			this.y = 0.25 * s;
			this.z = (m12 + m21) / s;

		} else {

			s = 2.0 * Math.sqrt(1.0 + m22 - m00 - m11);

			this.w = (m10 - m01) / s;
			this.x = (m02 + m20) / s;
			this.y = (m12 + m21) / s;
			this.z = 0.25 * s;

		}

		return this;

	};

	/**
		 * Sets the components of this quaternion based on unit vectors.
		 *
		 * @param {Vector3} vFrom - A unit vector. Assumed to be normalized.
		 * @param {Vector3} vTo - A unit vector. Assumed to be normalized.
		 * @return {Quaternion} This quaternion.
		 */

	Quaternion.prototype.setFromUnitVectors = function setFromUnitVectors (vFrom, vTo) {

		var r = vFrom.dot(vTo) + 1;

		if(r < 1e-6) {

			r = 0;

			if(Math.abs(vFrom.x) > Math.abs(vFrom.z)) {

				v$2.set(-vFrom.y, vFrom.x, 0);

			} else {

				v$2.set(0, -vFrom.z, vFrom.y);

			}

		} else {

			v$2.crossVectors(vFrom, vTo);

		}

		this.x = v$2.x;
		this.y = v$2.y;
		this.z = v$2.z;
		this.w = r;

		return this.normalize();

	};

	/**
		 * Inverts this quaternion.
		 *
		 * @return {Quaternion} This quaternion.
		 */

	Quaternion.prototype.invert = function invert () {

		return this.conjugate().normalize();

	};

	/**
		 * Conjugates this quaternion.
		 *
		 * @return {Quaternion} This quaternion.
		 */

	Quaternion.prototype.conjugate = function conjugate () {

		this.x *= -1;
		this.y *= -1;
		this.z *= -1;

		return this;

	};

	/**
		 * Calculates the squared length of this quaternion.
		 *
		 * @return {Number} The squared length.
		 */

	Quaternion.prototype.lengthSquared = function lengthSquared () {

		return this.x * this.x + this.y * this.y + this.z * this.z + this.w * this.w;

	};

	/**
		 * Calculates the length of this quaternion.
		 *
		 * @return {Number} The length.
		 */

	Quaternion.prototype.length = function length () {

		return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z + this.w * this.w);

	};

	/**
		 * Normalizes this quaternion.
		 *
		 * @return {Quaternion} This quaternion.
		 */

	Quaternion.prototype.normalize = function normalize () {

		var l = this.length();

		var invLength;

		if(l === 0) {

			this.x = 0;
			this.y = 0;
			this.z = 0;
			this.w = 1;

		} else {

			invLength = 1.0 / l;

			this.x = this.x * invLength;
			this.y = this.y * invLength;
			this.z = this.z * invLength;
			this.w = this.w * invLength;

		}

		return this;

	};

	/**
		 * Calculates the dot product with a given vector.
		 *
		 * @param {Vector4} v - A vector.
		 * @return {Number} The dot product.
		 */

	Quaternion.prototype.dot = function dot (v) {

		return this.x * v.x + this.y * v.y + this.z * v.z + this.w * v.w;

	};

	/**
		 * Multiplies the given quaternions and stores the result in this quaternion.
		 *
		 * For more details see:
		 *  http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/code/index.htm
		 *
		 * @param {Quaternion} q - A quaternion.
		 * @param {Quaternion} q - Another quaternion.
		 * @return {Quaternion} This quaternion.
		 */

	Quaternion.prototype.multiplyQuaternions = function multiplyQuaternions (a, b) {

		var qax = a.x, qay = a.y, qaz = a.z, qaw = a.w;
		var qbx = b.x, qby = b.y, qbz = b.z, qbw = b.w;

		this.x = qax * qbw + qaw * qbx + qay * qbz - qaz * qby;
		this.y = qay * qbw + qaw * qby + qaz * qbx - qax * qbz;
		this.z = qaz * qbw + qaw * qbz + qax * qby - qay * qbx;
		this.w = qaw * qbw - qax * qbx - qay * qby - qaz * qbz;

		return this;

	};

	/**
		 * Multiplies this quaternion with the given one and stores the result in
		 * this quaternion.
		 *
		 * @param {Quaternion} q - A quaternion.
		 * @return {Quaternion} This quaternion.
		 */

	Quaternion.prototype.multiply = function multiply (q) {

		return this.multiplyQuaternions(this, q);

	};

	/**
		 * Multiplies the given quaternion with this one and stores the result in
		 * this quaternion.
		 *
		 * @param {Quaternion} q - A quaternion.
		 * @return {Quaternion} This quaternion.
		 */

	Quaternion.prototype.premultiply = function premultiply (q) {

		return this.multiplyQuaternions(q, this);

	};

	/**
		 * Performs a spherical linear interpolation towards the given quaternion.
		 *
		 * For more details see:
		 *  http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/slerp/
		 *
		 * @param {Quaternion} q - A quaternion.
		 * @param {Number} t - The slerp factor.
		 * @return {Quaternion} This quaternion.
		 */

	Quaternion.prototype.slerp = function slerp (q, t) {

		var x = this.x, y = this.y, z = this.z, w = this.w;

		var cosHalfTheta, sinHalfTheta;
		var halfTheta, ratioA, ratioB;

		if(t === 1) {

			this.copy(q);

		} else if(t > 0) {

			cosHalfTheta = w * q.w + x * q.x + y * q.y + z * q.z;

			if(cosHalfTheta < 0) {

				this.w = -q.w;
				this.x = -q.x;
				this.y = -q.y;
				this.z = -q.z;

				cosHalfTheta = -cosHalfTheta;

			} else {

				this.copy(q);

			}

			if(cosHalfTheta >= 1.0) {

				this.w = w;
				this.x = x;
				this.y = y;
				this.z = z;

				return this;

			}

			sinHalfTheta = Math.sqrt(1.0 - cosHalfTheta * cosHalfTheta);

			if(Math.abs(sinHalfTheta) < 1e-3) {

				this.w = 0.5 * (w + this.w);
				this.x = 0.5 * (x + this.x);
				this.y = 0.5 * (y + this.y);
				this.z = 0.5 * (z + this.z);

				return this;

			}

			halfTheta = Math.atan2(sinHalfTheta, cosHalfTheta);
			ratioA = Math.sin((1.0 - t) * halfTheta) / sinHalfTheta;
			ratioB = Math.sin(t * halfTheta) / sinHalfTheta;

			this.w = (w * ratioA + this.w * ratioB);
			this.x = (x * ratioA + this.x * ratioB);
			this.y = (y * ratioA + this.y * ratioB);
			this.z = (z * ratioA + this.z * ratioB);

		}

		return this;

	};

	/**
		 * Checks if this quaternions equals the given one.
		 *
		 * @return {Boolean} Whether the quaternions are equal.
		 */

	Quaternion.prototype.equals = function equals (q) {

		return (q.x === this.x) && (q.y === this.y) && (q.z === this.z) && (q.w === this.w);

	};

	/**
		 * Performs a spherical linear interpolation.
		 *
		 * @param {Quaternion} qa - The base quaternion.
		 * @param {Quaternion} qb - The target quaternion.
		 * @param {Quaternion} qr - A quaternion to store the result in. If none is provided, a new one will be created.
		 * @param {Number} t - The slerp factor.
		 * @return {Quaternion} The resulting quaternion.
		 */

	Quaternion.slerp = function slerp (qa, qb, qr, t) {
			if ( qr === void 0 ) qr = new Quaternion();


		return qr.copy(qa).slerp(qb, t);

	};

	/**
		 * Performs an array-based spherical linear interpolation.
		 *
		 * @param {Number[]} dst - An array to store the result in.
		 * @param {Number} dstOffset - An offset into the destination array.
		 * @param {Number[]} src0 - An array that contains the base quaternion values.
		 * @param {Number} src0Offset - An offset into the base array.
		 * @param {Number[]} src1 - An array that contains the target quaternion values.
		 * @param {Number} src1Offset - An offset into the target array.
		 * @param {Number} t - The slerp factor.
		 */

	Quaternion.slerpFlat = function slerpFlat (dst, dstOffset, src0, srcOffset0, src1, srcOffset1, t) {

		var x1 = src1[srcOffset1];
		var y1 = src1[srcOffset1 + 1];
		var z1 = src1[srcOffset1 + 2];
		var w1 = src1[srcOffset1 + 3];

		var x0 = src0[srcOffset0];
		var y0 = src0[srcOffset0 + 1];
		var z0 = src0[srcOffset0 + 2];
		var w0 = src0[srcOffset0 + 3];

		var s, f;
		var sin, cos, sqrSin;
		var dir, len, tDir;

		if(w0 !== w1 || x0 !== x1 || y0 !== y1 || z0 !== z1) {

			s = 1.0 - t;
			cos = x0 * x1 + y0 * y1 + z0 * z1 + w0 * w1;

			dir = (cos >= 0) ? 1 : -1;
			sqrSin = 1.0 - cos * cos;

			// Skip the Slerp for tiny steps to avoid numeric problems.
			if(sqrSin > Number.EPSILON) {

				sin = Math.sqrt(sqrSin);
				len = Math.atan2(sin, cos * dir);

				s = Math.sin(s * len) / sin;
				t = Math.sin(t * len) / sin;

			}

			tDir = t * dir;

			x0 = x0 * s + x1 * tDir;
			y0 = y0 * s + y1 * tDir;
			z0 = z0 * s + z1 * tDir;
			w0 = w0 * s + w1 * tDir;

			// Normalize in case a lerp has just been performed.
			if(s === 1.0 - t) {

				f = 1.0 / Math.sqrt(x0 * x0 + y0 * y0 + z0 * z0 + w0 * w0);

				x0 *= f;
				y0 *= f;
				z0 *= f;
				w0 *= f;

			}

		}

		dst[dstOffset] = x0;
		dst[dstOffset + 1] = y0;
		dst[dstOffset + 2] = z0;
		dst[dstOffset + 3] = w0;

	};

	/**
	 * Clamps the given value.
	 *
	 * @private
	 * @param {Number} value - The value.
	 * @param {Number} min - The lower limit.
	 * @param {Number} max - The upper limit.
	 * @return {Number} The clamped value.
	 */

	function clamp(value, min, max) {

		return Math.max(Math.min(value, max), min);

	}

	/**
	 * A matrix.
	 *
	 * @type {Matrix3}
	 * @private
	 */

	var m = new Matrix3();

	/**
	 * A quaternion.
	 *
	 * @type {Quaternion}
	 * @private
	 */

	var q = new Quaternion();

	/**
	 * Euler angles.
	 */

	var Euler = function Euler(x, y, z) {
		if ( x === void 0 ) x = 0;
		if ( y === void 0 ) y = 0;
		if ( z === void 0 ) z = 0;


		/**
			 * The rotation around the X-axis.
			 *
			 * @type {Number}
			 */

		this.x = x;

		/**
			 * The rotation around the Y-axis.
			 *
			 * @type {Number}
			 */

		this.y = y;

		/**
			 * The rotation around the Z-axis.
			 *
			 * @type {Number}
			 */

		this.z = z;

		/**
			 * The rotation order.
			 *
			 * @type {RotationOrder}
			 * @default Euler.defaultOrder
			 */

		this.order = Euler.defaultOrder;

	};

	var staticAccessors = { defaultOrder: { configurable: true } };

	/**
		 * Sets the Euler angles and rotation order.
		 *
		 * @param {Number} x - The rotation around the X-axis.
		 * @param {Number} y - The rotation around the Y-axis.
		 * @param {Number} z - The rotation around the Z-axis.
		 * @param {Number} order - The rotation order.
		 * @return {Euler} This set of Euler angles.
		 */

	Euler.prototype.set = function set (x, y, z, order) {

		this.x = x;
		this.y = y;
		this.z = z;
		this.order = z;

		return this;

	};

	/**
		 * Copies the values of another set of Euler angles.
		 *
		 * @param {Euler} e - A set of Euler angles.
		 * @return {Euler} This set of Euler angles.
		 */

	Euler.prototype.copy = function copy (e) {

		this.x = e.x;
		this.y = e.y;
		this.z = e.z;
		this.order = e.order;

		return this;

	};

	/**
		 * Clones this set of Euler angles.
		 *
		 * @return {Euler} A clone of this set of Euler angles.
		 */

	Euler.prototype.clone = function clone () {

		return new this.constructor(this.x, this.y, this.z, this.order);

	};

	/**
		 * Copies angles and the rotation order from an array.
		 *
		 * @param {Number[]} array - An array.
		 * @param {Number} offset - An offset.
		 * @return {Euler} This set of Euler angles.
		 */

	Euler.prototype.fromArray = function fromArray (array, offset) {
			if ( offset === void 0 ) offset = 0;


		this.x = array[offset];
		this.y = array[offset + 1];
		this.z = array[offset + 2];
		this.order = array[offset + 3];

		return this;

	};

	/**
		 * Stores this set of Euler angles and the rotation order in an array.
		 *
		 * @param {Array} [array] - A target array.
		 * @param {Number} offset - An offset.
		 * @return {Number[]} The array.
		 */

	Euler.prototype.toArray = function toArray (array, offset) {
			if ( array === void 0 ) array = [];
			if ( offset === void 0 ) offset = 0;


		array[offset] = this.x;
		array[offset + 1] = this.y;
		array[offset + 2] = this.z;
		array[offset + 3] = this.order;

		return array;

	};

	/**
		 * Stores this set of Euler angles in a vector.
		 *
		 * @param {Vector3} [vector] - A target vector. If none is provided, a new one will be created.
		 * @return {Vector3} The vector.
		 */

	Euler.prototype.toVector3 = function toVector3 (vector) {
			if ( vector === void 0 ) vector = new Vector3();


		return vector.set(this.x, this.y, this.z);

	};

	/**
		 * Copies the rotation from a given matrix.
		 *
		 * @param {Matrix4} m - A rotation matrix. The upper 3x3 is assumed to be a pure rotation matrix (i.e. unscaled).
		 * @param {RotationOrder} [order] - An override rotation order.
		 * @return {Euler} This set of Euler angles.
		 */

	Euler.prototype.setFromRotationMatrix = function setFromRotationMatrix (m, order) {
			if ( order === void 0 ) order = this.order;


		var te = m.elements;
		var m00 = te[0], m01 = te[4], m02 = te[8];
		var m10 = te[1], m11 = te[5], m12 = te[9];
		var m20 = te[2], m21 = te[6], m22 = te[10];

		var THRESHOLD = 1.0 - 1e-5;

		switch(order) {

			case RotationOrder.XYZ: {

				this.y = Math.asin(clamp(m02, -1, 1));

				if(Math.abs(m02) < THRESHOLD) {

					this.x = Math.atan2(-m12, m22);
					this.z = Math.atan2(-m01, m00);

				} else {

					this.x = Math.atan2(m21, m11);
					this.z = 0;

				}

				break;

			}

			case RotationOrder.YXZ: {

				this.x = Math.asin(-clamp(m12, -1, 1));

				if(Math.abs(m12) < THRESHOLD) {

					this.y = Math.atan2(m02, m22);
					this.z = Math.atan2(m10, m11);

				} else {

					this.y = Math.atan2(-m20, m00);
					this.z = 0;

				}

				break;

			}

			case RotationOrder.ZXY: {

				this.x = Math.asin(clamp(m21, -1, 1));

				if(Math.abs(m21) < THRESHOLD) {

					this.y = Math.atan2(-m20, m22);
					this.z = Math.atan2(-m01, m11);

				} else {

					this.y = 0;
					this.z = Math.atan2(m10, m00);

				}

				break;

			}

			case RotationOrder.ZYX: {

				this.y = Math.asin(-clamp(m20, -1, 1));

				if(Math.abs(m20) < THRESHOLD) {

					this.x = Math.atan2(m21, m22);
					this.z = Math.atan2(m10, m00);

				} else {

					this.x = 0;
					this.z = Math.atan2(-m01, m11);

				}

				break;

			}

			case RotationOrder.YZX: {

				this.z = Math.asin(clamp(m10, -1, 1));

				if(Math.abs(m10) < THRESHOLD) {

					this.x = Math.atan2(-m12, m11);
					this.y = Math.atan2(-m20, m00);

				} else {

					this.x = 0;
					this.y = Math.atan2(m02, m22);

				}

				break;

			}

			case RotationOrder.XZY: {

				this.z = Math.asin(-clamp(m01, -1, 1));

				if(Math.abs(m01) < THRESHOLD) {

					this.x = Math.atan2(m21, m11);
					this.y = Math.atan2(m02, m00);

				} else {

					this.x = Math.atan2(-m12, m22);
					this.y = 0;

				}

				break;

			}

		}

		this.order = order;

		return this;

	};

	/**
		 * Copies the rotation from a given quaternion.
		 *
		 * @param {Matrix4} q - A quaternion.
		 * @param {RotationOrder} [order] - An override rotation order.
		 * @return {Euler} This set of Euler angles.
		 */

	Euler.prototype.setFromQuaternion = function setFromQuaternion (q, order) {

		m.makeRotationFromQuaternion(q);

		return this.setFromRotationMatrix(m, order);

	};

	/**
		 * Copies the rotation from a given vector.
		 *
		 * @param {Matrix4} v - A vector.
		 * @param {RotationOrder} [order] - A rotation order.
		 * @return {Euler} This set of Euler angles.
		 */

	Euler.prototype.setFromVector3 = function setFromVector3 (v, order) {
			if ( order === void 0 ) order = this.order;


		return this.set(v.x, v.y, v.z, order);

	};

	/**
		 * Reorder the rotation angles.
		 *
		 * WARNING: this operation discards revolution information!
		 *
		 * @param {RotationOrder} newOrder - The new rotation order.
		 * @return {Euler} This set of Euler angles.
		 */

	Euler.prototype.reorder = function reorder (newOrder) {

		q.setFromEuler(this);

		return this.setFromQuaternion(q, newOrder);

	};

	/**
		 * Checks if this set of Euler angles equals the given one.
		 *
		 * @param {Euler} e - Euler angles.
		 * @return {Boolean} Whether this set of Euler angles equals the given one.
		 */

	Euler.prototype.equals = function equals (e) {

		return (e.x === this.x && e.y === this.y && e.z === this.z && e.order === this.order);

	};

	/**
		 * The default rotation order.
		 *
		 * @type {RotationOrder}
		 * @final
		 */

	staticAccessors.defaultOrder.get = function () { return RotationOrder.XYZ; };

	Object.defineProperties( Euler, staticAccessors );

	/**
	 * A vector.
	 *
	 * @type {Vector3}
	 * @private
	 */

	var a = new Vector3();

	/**
	 * A vector.
	 *
	 * @type {Vector3}
	 * @private
	 */

	var b = new Vector3();

	/**
	 * A line.
	 */

	var Line3 = function Line3(start, end) {
		if ( start === void 0 ) start = new Vector3();
		if ( end === void 0 ) end = new Vector3();


		/**
			 * The starting point.
			 *
			 * @type {Vector3}
			 */

		this.start = start;

		/**
			 * The ending point.
			 *
			 * @type {Vector3}
			 */

		this.end = end;

	};

	/**
		 * Sets the starting and ending point of this line.
		 *
		 * @param {Vector3} start - The starting point.
		 * @param {Vector3} end - The ending point.
		 * @return {Line3} This line.
		 */

	Line3.prototype.set = function set (start, end) {

		this.start.copy(start);
		this.end.copy(end);

		return this;

	};

	/**
		 * Copies the values of the given line.
		 *
		 * @param {Line3} l - A line.
		 * @return {Line3} This line.
		 */

	Line3.prototype.copy = function copy (l) {

		this.start.copy(l.start);
		this.end.copy(l.end);

		return this;

	};

	/**
		 * Clones this line.
		 *
		 * @return {Line3} The cloned line.
		 */

	Line3.prototype.clone = function clone () {

		return new this.constructor().copy(this);

	};

	/**
		 * Calculates the center of this line.
		 *
		 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
		 * @return {Vector3} The center of this line.
		 */

	Line3.prototype.getCenter = function getCenter (target) {
			if ( target === void 0 ) target = new Vector3();


		return target.addVectors(this.start, this.end).multiplyScalar(0.5);

	};

	/**
		 * Calculates the delta vector of this line.
		 *
		 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
		 * @return {Vector3} The delta vector of this line.
		 */

	Line3.prototype.delta = function delta (target) {
			if ( target === void 0 ) target = new Vector3();


		return target.subVectors(this.end, this.start);

	};

	/**
		 * Calculates the squared length of this line.
		 *
		 * @return {Vector3} The squared length.
		 */

	Line3.prototype.lengthSquared = function lengthSquared () {

		return this.start.distanceToSquared(this.end);

	};

	/**
		 * Calculates the length of this line.
		 *
		 * @return {Vector3} The length.
		 */

	Line3.prototype.length = function length () {

		return this.start.distanceTo(this.end);

	};

	/**
		 * Adjusts the lin to point in the given direction.
		 *
		 * @param {Vector3} d - The direction.
		 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
		 * @return {Vector3} The length.
		 */

	Line3.prototype.at = function at (d, target) {

		return this.delta(target).multiplyScalar(d).add(this.start);

	};

	/**
		 * Returns a point parameter based on the closest point as projected on the line segement.
		 *
		 * @private
		 * @param {Vector3} p - A point.
		 * @param {Boolean} clampToLine - Whether the point should be clamped to the line.
		 * @return {Vector3} The parameter.
		 */

	Line3.prototype.closestPointToPointParameter = function closestPointToPointParameter (p, clampToLine) {

		a.subVectors(p, this.start);
		b.subVectors(this.end, this.start);

		var bb = b.dot(b);
		var ba = b.dot(a);

		var t = clampToLine ? Math.min(Math.max(ba / bb, 0), 1) : ba / bb;

		return t;

	};

	/**
		 * Returns the closets point on the line.
		 *
		 * @param {Vector3} p - A point.
		 * @param {Boolean} [clampToLine=false] - Whether the point should be clamped to the line.
		 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
		 * @return {Vector3} The parameter.
		 */

	Line3.prototype.closestPointToPoint = function closestPointToPoint (p, clampToLine, target) {
			if ( clampToLine === void 0 ) clampToLine = false;
			if ( target === void 0 ) target = new Vector3();


		var t = this.closestPointToPointParameter(p, clampToLine);

		return this.delta(target).multiplyScalar(t).add(this.start);

	};

	/**
		 * Checks if this line equals the given one.
		 *
		 * @param {Line3} l - A line.
		 * @return {Boolean} Whether the lines are equal.
		 */

	Line3.prototype.equals = function equals (l) {

		return l.start.equals(this.start) && l.end.equals(this.end);

	};

	/**
	 * A vector.
	 *
	 * @type {Vector3}
	 * @private
	 */

	var a$1 = new Vector3();

	/**
	 * A vector.
	 *
	 * @type {Vector3}
	 * @private
	 */

	var b$1 = new Vector3();

	/**
	 * A vector.
	 *
	 * @type {Vector3}
	 * @private
	 */

	var c$1 = new Vector3();

	/**
	 * A 4x4 matrix.
	 */

	var Matrix4 = function Matrix4() {

		/**
			 * The matrix elements.
			 *
			 * @type {Float32Array}
			 */

		this.elements = new Float32Array([

			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1

		]);

	};

	/**
		 * Sets the values of this matrix.
		 *
		 * @param {Number} n00 - The value of the first row, first column.
		 * @param {Number} n01 - The value of the first row, second column.
		 * @param {Number} n02 - The value of the first row, third column.
		 * @param {Number} n03 - The value of the first row, fourth column.
		 * @param {Number} n10 - The value of the second row, first column.
		 * @param {Number} n11 - The value of the second row, second column.
		 * @param {Number} n12 - The value of the second row, third column.
		 * @param {Number} n13 - The value of the second row, fourth column.
		 * @param {Number} n20 - The value of the third row, first column.
		 * @param {Number} n21 - The value of the third row, second column.
		 * @param {Number} n22 - The value of the third row, third column.
		 * @param {Number} n23 - The value of the third row, fourth column.
		 * @param {Number} n30 - The value of the fourth row, first column.
		 * @param {Number} n31 - The value of the fourth row, second column.
		 * @param {Number} n32 - The value of the fourth row, third column.
		 * @param {Number} n33 - The value of the fourth row, fourth column.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.set = function set (n00, n01, n02, n03, n10, n11, n12, n13, n20, n21, n22, n23, n30, n31, n32, n33) {

		var te = this.elements;

		te[0] = n00; te[4] = n01; te[8] = n02; te[12] = n03;
		te[1] = n10; te[5] = n11; te[9] = n12; te[13] = n13;
		te[2] = n20; te[6] = n21; te[10] = n22; te[14] = n23;
		te[3] = n30; te[7] = n31; te[11] = n32; te[15] = n33;

		return this;

	};

	/**
		 * Sets this matrix to the identity matrix.
		 *
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.identity = function identity () {

		this.set(

			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1

		);

		return this;

	};

	/**
		 * Copies the values of a given matrix.
		 *
		 * @param {Matrix4} matrix - A matrix.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.copy = function copy (matrix) {

		var me = matrix.elements;
		var te = this.elements;

		te[0] = me[0]; te[1] = me[1]; te[2] = me[2]; te[3] = me[3];
		te[4] = me[4]; te[5] = me[5]; te[6] = me[6]; te[7] = me[7];
		te[8] = me[8]; te[9] = me[9]; te[10] = me[10]; te[11] = me[11];
		te[12] = me[12]; te[13] = me[13]; te[14] = me[14]; te[15] = me[15];

		return this;

	};

	/**
		 * Clones this matrix.
		 *
		 * @return {Matrix4} A clone of this matrix.
		 */

	Matrix4.prototype.clone = function clone () {

		return new this.constructor().fromArray(this.elements);

	};

	/**
		 * Copies the values of a given array.
		 *
		 * @param {Number[]} array - An array.
		 * @param {Number} [offset=0] - An offset into the array.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.fromArray = function fromArray (array, offset) {
			if ( offset === void 0 ) offset = 0;


		var te = this.elements;

		var i;

		for(i = 0; i < 16; ++i) {

			te[i] = array[i + offset];

		}

		return this;

	};

	/**
		 * Stores this matrix in an array.
		 *
		 * @param {Number[]} [array] - A target array.
		 * @param {Number} [offset=0] - An offset into the array.
		 * @return {Number[]} The array.
		 */

	Matrix4.prototype.toArray = function toArray (array, offset) {
			if ( array === void 0 ) array = [];
			if ( offset === void 0 ) offset = 0;


		var te = this.elements;

		var i;

		for(i = 0; i < 16; ++i) {

			array[i + offset] = te[i];

		}

		return array;

	};

	/**
		 * Returns the largest scale.
		 *
		 * @param {Matrix4} matrix - A matrix.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.getMaxScaleOnAxis = function getMaxScaleOnAxis () {

		var te = this.elements;

		var scaleXSq = te[0] * te[0] + te[1] * te[1] + te[2] * te[2];
		var scaleYSq = te[4] * te[4] + te[5] * te[5] + te[6] * te[6];
		var scaleZSq = te[8] * te[8] + te[9] * te[9] + te[10] * te[10];

		return Math.sqrt(Math.max(scaleXSq, scaleYSq, scaleZSq));

	};

	/**
		 * Copies the position values of a given matrix.
		 *
		 * @param {Matrix4} matrix - A matrix.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.copyPosition = function copyPosition (matrix) {

		var te = this.elements;
		var me = matrix.elements;

		te[12] = me[12];
		te[13] = me[13];
		te[14] = me[14];

		return this;

	};

	/**
		 * Sets the position values of this matrix.
		 *
		 * @param {Vector3} p - A position.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.setPosition = function setPosition (p) {

		var te = this.elements;

		te[12] = p.x;
		te[13] = p.y;
		te[14] = p.z;

		return this;

	};

	/**
		 * Extracts the basis from this matrix.
		 *
		 * @param {Vector3} xAxis - A vector to store the X-axis column in.
		 * @param {Vector3} yAxis - A vector to store the Y-axis column in.
		 * @param {Vector3} zAxis - A vector to store the Z-axis column in.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.extractBasis = function extractBasis (xAxis, yAxis, zAxis) {

		xAxis.setFromMatrixColumn(this, 0);
		yAxis.setFromMatrixColumn(this, 1);
		zAxis.setFromMatrixColumn(this, 2);

		return this;

	};

	/**
		 * Sets the basis of this matrix.
		 *
		 * @param {Vector3} xAxis - The X-axis.
		 * @param {Vector3} yAxis - The Y-axis.
		 * @param {Vector3} zAxis - The Z-axis.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.makeBasis = function makeBasis (xAxis, yAxis, zAxis) {

		this.set(

			xAxis.x, yAxis.x, zAxis.x, 0,
			xAxis.y, yAxis.y, zAxis.y, 0,
			xAxis.z, yAxis.z, zAxis.z, 0,
			0, 0, 0, 1

		);

		return this;

	};

	/**
		 * Extracts the rotation from a given matrix.
		 *
		 * @param {Matrix4} m - A matrix.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.extractRotation = function extractRotation (m) {

		var te = this.elements;
		var me = m.elements;

		var scaleX = 1.0 / a$1.setFromMatrixColumn(m, 0).length();
		var scaleY = 1.0 / a$1.setFromMatrixColumn(m, 1).length();
		var scaleZ = 1.0 / a$1.setFromMatrixColumn(m, 2).length();

		te[0] = me[0] * scaleX;
		te[1] = me[1] * scaleX;
		te[2] = me[2] * scaleX;

		te[4] = me[4] * scaleY;
		te[5] = me[5] * scaleY;
		te[6] = me[6] * scaleY;

		te[8] = me[8] * scaleZ;
		te[9] = me[9] * scaleZ;
		te[10] = me[10] * scaleZ;

		return this;

	};

	/**
		 * Sets the matrix rotation based on the given Euler angles.
		 *
		 * @param {Euler} euler - The euler angles.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.makeRotationFromEuler = function makeRotationFromEuler (euler) {

		var te = this.elements;

		var x = euler.x;
		var y = euler.y;
		var z = euler.z;

		var a = Math.cos(x), b = Math.sin(x);
		var c = Math.cos(y), d = Math.sin(y);
		var e = Math.cos(z), f = Math.sin(z);

		var ae, af, be, bf;
		var ce, cf, de, df;
		var ac, ad, bc, bd;

		switch(euler.order) {

			case RotationOrder.XYZ: {

				ae = a * e, af = a * f, be = b * e, bf = b * f;

				te[0] = c * e;
				te[4] = -c * f;
				te[8] = d;

				te[1] = af + be * d;
				te[5] = ae - bf * d;
				te[9] = -b * c;

				te[2] = bf - ae * d;
				te[6] = be + af * d;
				te[10] = a * c;

				break;

			}

			case RotationOrder.YXZ: {

				ce = c * e, cf = c * f, de = d * e, df = d * f;

				te[0] = ce + df * b;
				te[4] = de * b - cf;
				te[8] = a * d;

				te[1] = a * f;
				te[5] = a * e;
				te[9] = -b;

				te[2] = cf * b - de;
				te[6] = df + ce * b;
				te[10] = a * c;

				break;

			}

			case RotationOrder.ZXY: {

				ce = c * e, cf = c * f, de = d * e, df = d * f;

				te[0] = ce - df * b;
				te[4] = -a * f;
				te[8] = de + cf * b;

				te[1] = cf + de * b;
				te[5] = a * e;
				te[9] = df - ce * b;

				te[2] = -a * d;
				te[6] = b;
				te[10] = a * c;

				break;

			}

			case RotationOrder.ZYX: {

				ae = a * e, af = a * f, be = b * e, bf = b * f;

				te[0] = c * e;
				te[4] = be * d - af;
				te[8] = ae * d + bf;

				te[1] = c * f;
				te[5] = bf * d + ae;
				te[9] = af * d - be;

				te[2] = -d;
				te[6] = b * c;
				te[10] = a * c;

				break;

			}

			case RotationOrder.YZX: {

				ac = a * c, ad = a * d, bc = b * c, bd = b * d;

				te[0] = c * e;
				te[4] = bd - ac * f;
				te[8] = bc * f + ad;

				te[1] = f;
				te[5] = a * e;
				te[9] = -b * e;

				te[2] = -d * e;
				te[6] = ad * f + bc;
				te[10] = ac - bd * f;

				break;

			}

			case RotationOrder.XZY: {

				ac = a * c, ad = a * d, bc = b * c, bd = b * d;

				te[0] = c * e;
				te[4] = -f;
				te[8] = d * e;

				te[1] = ac * f + bd;
				te[5] = a * e;
				te[9] = ad * f - bc;

				te[2] = bc * f - ad;
				te[6] = b * e;
				te[10] = bd * f + ac;

				break;

			}

		}

		// Last column.
		te[3] = 0;
		te[7] = 0;
		te[11] = 0;

		// Bottom row.
		te[12] = 0;
		te[13] = 0;
		te[14] = 0;
		te[15] = 1;

		return this;

	};

	/**
		 * Sets the matrix rotation based on the given quaternion.
		 *
		 * @param {Quaternion} q - The quaternion.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.makeRotationFromQuaternion = function makeRotationFromQuaternion (q) {

		var te = this.elements;

		var x = q.x, y = q.y, z = q.z, w = q.w;
		var x2 = x + x, y2 = y + y, z2 = z + z;
		var xx = x * x2, xy = x * y2, xz = x * z2;
		var yy = y * y2, yz = y * z2, zz = z * z2;
		var wx = w * x2, wy = w * y2, wz = w * z2;

		te[0] = 1 - (yy + zz);
		te[4] = xy - wz;
		te[8] = xz + wy;

		te[1] = xy + wz;
		te[5] = 1 - (xx + zz);
		te[9] = yz - wx;

		te[2] = xz - wy;
		te[6] = yz + wx;
		te[10] = 1 - (xx + yy);

		// Last column.
		te[3] = 0;
		te[7] = 0;
		te[11] = 0;

		// Bottom row.
		te[12] = 0;
		te[13] = 0;
		te[14] = 0;
		te[15] = 1;

		return this;

	};

	/**
		 * Creates a rotation that looks at the given target.
		 *
		 * @param {Vector3} eye - The position of the eye.
		 * @param {Vector3} target - The target to look at.
		 * @param {Vector3} up - The up vector.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.lookAt = function lookAt (eye, target, up) {

		var te = this.elements;
		var x = a$1, y = b$1, z = c$1;

		z.subVectors(eye, target);

		if(z.lengthSquared() === 0) {

			// Eye and target are at the same position.
			z.z = 1;

		}

		z.normalize();
		x.crossVectors(up, z);

		if(x.lengthSquared() === 0) {

			// Up and z are parallel.
			if(Math.abs(up.z) === 1) {

				z.x += 1e-4;

			} else {

				z.z += 1e-4;

			}

			z.normalize();
			x.crossVectors(up, z);

		}

		x.normalize();
		y.crossVectors(z, x);

		te[0] = x.x; te[4] = y.x; te[8] = z.x;
		te[1] = x.y; te[5] = y.y; te[9] = z.y;
		te[2] = x.z; te[6] = y.z; te[10] = z.z;

		return this;

	};

	/**
		 * Sets this matrix to the product of the given matrices.
		 *
		 * @param {Matrix4} a - A matrix.
		 * @param {Matrix4} b - A matrix.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.multiplyMatrices = function multiplyMatrices (a, b) {

		var te = this.elements;
		var ae = a.elements;
		var be = b.elements;

		var a00 = ae[0], a01 = ae[4], a02 = ae[8], a03 = ae[12];
		var a10 = ae[1], a11 = ae[5], a12 = ae[9], a13 = ae[13];
		var a20 = ae[2], a21 = ae[6], a22 = ae[10], a23 = ae[14];
		var a30 = ae[3], a31 = ae[7], a32 = ae[11], a33 = ae[15];

		var b00 = be[0], b01 = be[4], b02 = be[8], b03 = be[12];
		var b10 = be[1], b11 = be[5], b12 = be[9], b13 = be[13];
		var b20 = be[2], b21 = be[6], b22 = be[10], b23 = be[14];
		var b30 = be[3], b31 = be[7], b32 = be[11], b33 = be[15];

		te[0] = a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30;
		te[4] = a00 * b01 + a01 * b11 + a02 * b21 + a03 * b31;
		te[8] = a00 * b02 + a01 * b12 + a02 * b22 + a03 * b32;
		te[12] = a00 * b03 + a01 * b13 + a02 * b23 + a03 * b33;

		te[1] = a10 * b00 + a11 * b10 + a12 * b20 + a13 * b30;
		te[5] = a10 * b01 + a11 * b11 + a12 * b21 + a13 * b31;
		te[9] = a10 * b02 + a11 * b12 + a12 * b22 + a13 * b32;
		te[13] = a10 * b03 + a11 * b13 + a12 * b23 + a13 * b33;

		te[2] = a20 * b00 + a21 * b10 + a22 * b20 + a23 * b30;
		te[6] = a20 * b01 + a21 * b11 + a22 * b21 + a23 * b31;
		te[10] = a20 * b02 + a21 * b12 + a22 * b22 + a23 * b32;
		te[14] = a20 * b03 + a21 * b13 + a22 * b23 + a23 * b33;

		te[3] = a30 * b00 + a31 * b10 + a32 * b20 + a33 * b30;
		te[7] = a30 * b01 + a31 * b11 + a32 * b21 + a33 * b31;
		te[11] = a30 * b02 + a31 * b12 + a32 * b22 + a33 * b32;
		te[15] = a30 * b03 + a31 * b13 + a32 * b23 + a33 * b33;

		return this;

	};

	/**
		 * Multiplies this matrix with a given one.
		 *
		 * @param {Matrix4} m - A matrix.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.multiply = function multiply (m) {

		return this.multiplyMatrices(this, m);

	};

	/**
		 * Multiplies a given matrix with this one.
		 *
		 * @param {Matrix4} m - A matrix.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.premultiply = function premultiply (m) {

		return this.multiplyMatrices(m, this);

	};

	/**
		 * Multiplies this matrix with a given scalar.
		 *
		 * @param {Number} m - A scalar.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.multiplyScalar = function multiplyScalar (s) {

		var te = this.elements;

		te[0] *= s; te[4] *= s; te[8] *= s; te[12] *= s;
		te[1] *= s; te[5] *= s; te[9] *= s; te[13] *= s;
		te[2] *= s; te[6] *= s; te[10] *= s; te[14] *= s;
		te[3] *= s; te[7] *= s; te[11] *= s; te[15] *= s;

		return this;

	};

	/**
		 * Calculates the determinant of this matrix.
		 *
		 * For more details see:
		 *  http://www.euclideanspace.com/maths/algebra/matrix/functions/inverse/fourD/index.htm
		 *
		 * @return {Number} The determinant.
		 */

	Matrix4.prototype.determinant = function determinant () {

		var te = this.elements;

		var n00 = te[0], n01 = te[4], n02 = te[8], n03 = te[12];
		var n10 = te[1], n11 = te[5], n12 = te[9], n13 = te[13];
		var n20 = te[2], n21 = te[6], n22 = te[10], n23 = te[14];
		var n30 = te[3], n31 = te[7], n32 = te[11], n33 = te[15];

		var n00n11 = n00 * n11, n00n12 = n00 * n12, n00n13 = n00 * n13;
		var n01n10 = n01 * n10, n01n12 = n01 * n12, n01n13 = n01 * n13;
		var n02n10 = n02 * n10, n02n11 = n02 * n11, n02n13 = n02 * n13;
		var n03n10 = n03 * n10, n03n11 = n03 * n11, n03n12 = n03 * n12;

		return (

			n30 * (
				n03n12 * n21 -
				n02n13 * n21 -
				n03n11 * n22 +
				n01n13 * n22 +
				n02n11 * n23 -
				n01n12 * n23
			) +

			n31 * (
				n00n12 * n23 -
				n00n13 * n22 +
				n03n10 * n22 -
				n02n10 * n23 +
				n02n13 * n20 -
				n03n12 * n20
			) +

			n32 * (
				n00n13 * n21 -
				n00n11 * n23 -
				n03n10 * n21 +
				n01n10 * n23 +
				n03n11 * n20 -
				n01n13 * n20
			) +

			n33 * (
				-n02n11 * n20 -
				n00n12 * n21 +
				n00n11 * n22 +
				n02n10 * n21 -
				n01n10 * n22 +
				n01n12 * n20
			)

		);

	};

	/**
		 * Inverts the given matrix and stores the result in this matrix.
		 *
		 * For details see:
		 *  http://www.euclideanspace.com/maths/algebra/matrix/functions/inverse/fourD/index.htm
		 *
		 * @param {Matrix4} matrix - The matrix that should be inverted.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.getInverse = function getInverse (matrix) {

		var te = this.elements;
		var me = matrix.elements;

		var n00 = me[0], n10 = me[1], n20 = me[2], n30 = me[3];
		var n01 = me[4], n11 = me[5], n21 = me[6], n31 = me[7];
		var n02 = me[8], n12 = me[9], n22 = me[10], n32 = me[11];
		var n03 = me[12], n13 = me[13], n23 = me[14], n33 = me[15];

		var t00 = n12 * n23 * n31 - n13 * n22 * n31 + n13 * n21 * n32 - n11 * n23 * n32 - n12 * n21 * n33 + n11 * n22 * n33;
		var t01 = n03 * n22 * n31 - n02 * n23 * n31 - n03 * n21 * n32 + n01 * n23 * n32 + n02 * n21 * n33 - n01 * n22 * n33;
		var t02 = n02 * n13 * n31 - n03 * n12 * n31 + n03 * n11 * n32 - n01 * n13 * n32 - n02 * n11 * n33 + n01 * n12 * n33;
		var t03 = n03 * n12 * n21 - n02 * n13 * n21 - n03 * n11 * n22 + n01 * n13 * n22 + n02 * n11 * n23 - n01 * n12 * n23;

		var det = n00 * t00 + n10 * t01 + n20 * t02 + n30 * t03;

		var invDet;

		if(det !== 0) {

			invDet = 1.0 / det;

			te[0] = t00 * invDet;
			te[1] = (n13 * n22 * n30 - n12 * n23 * n30 - n13 * n20 * n32 + n10 * n23 * n32 + n12 * n20 * n33 - n10 * n22 * n33) * invDet;
			te[2] = (n11 * n23 * n30 - n13 * n21 * n30 + n13 * n20 * n31 - n10 * n23 * n31 - n11 * n20 * n33 + n10 * n21 * n33) * invDet;
			te[3] = (n12 * n21 * n30 - n11 * n22 * n30 - n12 * n20 * n31 + n10 * n22 * n31 + n11 * n20 * n32 - n10 * n21 * n32) * invDet;

			te[4] = t01 * invDet;
			te[5] = (n02 * n23 * n30 - n03 * n22 * n30 + n03 * n20 * n32 - n00 * n23 * n32 - n02 * n20 * n33 + n00 * n22 * n33) * invDet;
			te[6] = (n03 * n21 * n30 - n01 * n23 * n30 - n03 * n20 * n31 + n00 * n23 * n31 + n01 * n20 * n33 - n00 * n21 * n33) * invDet;
			te[7] = (n01 * n22 * n30 - n02 * n21 * n30 + n02 * n20 * n31 - n00 * n22 * n31 - n01 * n20 * n32 + n00 * n21 * n32) * invDet;

			te[8] = t02 * invDet;
			te[9] = (n03 * n12 * n30 - n02 * n13 * n30 - n03 * n10 * n32 + n00 * n13 * n32 + n02 * n10 * n33 - n00 * n12 * n33) * invDet;
			te[10] = (n01 * n13 * n30 - n03 * n11 * n30 + n03 * n10 * n31 - n00 * n13 * n31 - n01 * n10 * n33 + n00 * n11 * n33) * invDet;
			te[11] = (n02 * n11 * n30 - n01 * n12 * n30 - n02 * n10 * n31 + n00 * n12 * n31 + n01 * n10 * n32 - n00 * n11 * n32) * invDet;

			te[12] = t03 * invDet;
			te[13] = (n02 * n13 * n20 - n03 * n12 * n20 + n03 * n10 * n22 - n00 * n13 * n22 - n02 * n10 * n23 + n00 * n12 * n23) * invDet;
			te[14] = (n03 * n11 * n20 - n01 * n13 * n20 - n03 * n10 * n21 + n00 * n13 * n21 + n01 * n10 * n23 - n00 * n11 * n23) * invDet;
			te[15] = (n01 * n12 * n20 - n02 * n11 * n20 + n02 * n10 * n21 - n00 * n12 * n21 - n01 * n10 * n22 + n00 * n11 * n22) * invDet;

		} else {

			console.error("Can't invert matrix, determinant is zero", matrix);

			this.identity();

		}

		return this;

	};

	/**
		 * Transposes this matrix.
		 *
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.transpose = function transpose () {

		var te = this.elements;

		var t;

		t = te[1]; te[1] = te[4]; te[4] = t;
		t = te[2]; te[2] = te[8]; te[8] = t;
		t = te[6]; te[6] = te[9]; te[9] = t;

		t = te[3]; te[3] = te[12]; te[12] = t;
		t = te[7]; te[7] = te[13]; te[13] = t;
		t = te[11]; te[11] = te[14]; te[14] = t;

		return this;

	};

	/**
		 * Scales this matrix.
		 *
		 * @param {Number} sx - The X scale.
		 * @param {Number} sy - The Y scale.
		 * @param {Number} sy - The Z scale.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.scale = function scale (sx, sy, sz) {

		var te = this.elements;

		te[0] *= sx; te[4] *= sy; te[8] *= sz;
		te[1] *= sx; te[5] *= sy; te[9] *= sz;
		te[2] *= sx; te[6] *= sy; te[10] *= sz;
		te[3] *= sx; te[7] *= sy; te[11] *= sz;

		return this;

	};

	/**
		 * Makes this matrix a scale matrix.
		 *
		 * @param {Number} x - The X scale.
		 * @param {Number} y - The Y scale.
		 * @param {Number} z - The Z scale.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.makeScale = function makeScale (x, y, z) {

		this.set(

			x, 0, 0, 0,
			0, y, 0, 0,
			0, 0, z, 0,
			0, 0, 0, 1

		);

		return this;

	};

	/**
		 * Makes this matrix a translation matrix.
		 *
		 * @param {Number} x - The X offset.
		 * @param {Number} y - The Y offset.
		 * @param {Number} z - The Z offset.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.makeTranslation = function makeTranslation (x, y, z) {

		this.set(

			1, 0, 0, x,
			0, 1, 0, y,
			0, 0, 1, z,
			0, 0, 0, 1

		);

		return this;

	};

	/**
		 * Makes this matrix a rotation matrix.
		 *
		 * @param {Number} theta - The angle in radians.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.makeRotationX = function makeRotationX (theta) {

		var c = Math.cos(theta), s = Math.sin(theta);

		this.set(

			1, 0, 0, 0,
			0, c, -s, 0,
			0, s, c, 0,
			0, 0, 0, 1

		);

		return this;

	};

	/**
		 * Makes this matrix a rotation matrix with respect to the Y-axis.
		 *
		 * @param {Number} theta - The angle in radians.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.makeRotationY = function makeRotationY (theta) {

		var c = Math.cos(theta), s = Math.sin(theta);

		this.set(

			c, 0, s, 0,
			0, 1, 0, 0,
			-s, 0, c, 0,
			0, 0, 0, 1

		);

		return this;

	};

	/**
		 * Makes this matrix a rotation matrix with respect to the Z-axis.
		 *
		 * @param {Number} theta - The angle in radians.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.makeRotationZ = function makeRotationZ (theta) {

		var c = Math.cos(theta), s = Math.sin(theta);

		this.set(

			c, -s, 0, 0,
			s, c, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1

		);

		return this;

	};

	/**
		 * Makes this matrix a translation matrix with respect to a specific axis.
		 *
		 * For mor einformation see:
		 *  http://www.gamedev.net/reference/articles/article1199.asp
		 *
		 * @param {Vector3} axis - The axis. Assumed to be normalized.
		 * @param {Number} angle - The angle in radians.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.makeRotationAxis = function makeRotationAxis (axis, angle) {

		var c = Math.cos(angle);
		var s = Math.sin(angle);

		var t = 1.0 - c;

		var x = axis.x, y = axis.y, z = axis.z;
		var tx = t * x, ty = t * y;

		this.set(

			tx * x + c, tx * y - s * z, tx * z + s * y, 0,
			tx * y + s * z, ty * y + c, ty * z - s * x, 0,
			tx * z - s * y, ty * z + s * x, t * z * z + c, 0,
			0, 0, 0, 1

		);

		return this;

	};

	/**
		 * Makes this matrix a shear matrix.
		 *
		 * @param {Number} x - The X shear value.
		 * @param {Number} y - The Y shear value.
		 * @param {Number} z - The Z shear value.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.makeShear = function makeShear (x, y, z) {

		this.set(

			1, y, z, 0,
			x, 1, z, 0,
			x, y, 1, 0,
			0, 0, 0, 1

		);

		return this;

	};

	/**
		 * Sets this matrix based on the given position, rotation and scale.
		 *
		 * @param {Vector3} position - The position.
		 * @param {Quaternion} quaternion - The rotation.
		 * @param {Vector3} scale - The scale.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.compose = function compose (position, quaternion, scale) {

		this.makeRotationFromQuaternion(quaternion);
		this.scale(scale.x, scale.y, scale.z);
		this.setPosition(position);

		return this;

	};

	/**
		 * Decomposes this matrix into a position, rotation and scale vector.
		 *
		 * @param {Vector3} position - The target position.
		 * @param {Quaternion} quaternion - The target rotation.
		 * @param {Vector3} scale - The target scale.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.decompose = function decompose (position, quaternion, scale) {

		var te = this.elements;

		var n00 = te[0], n10 = te[1], n20 = te[2];
		var n01 = te[4], n11 = te[5], n21 = te[6];
		var n02 = te[8], n12 = te[9], n22 = te[10];

		var det = this.determinant();

		// If the determinant is negative, one scale must be inverted.
		var sx = a$1.set(n00, n10, n20).length() * ((det < 0) ? -1 : 1);
		var sy = a$1.set(n01, n11, n21).length();
		var sz = a$1.set(n02, n12, n22).length();

		var invSX = 1.0 / sx;
		var invSY = 1.0 / sy;
		var invSZ = 1.0 / sz;

		// Export the position.
		position.x = te[12];
		position.y = te[13];
		position.z = te[14];

		// Scale the rotation part.
		te[0] *= invSX; te[1] *= invSX; te[2] *= invSX;
		te[4] *= invSY; te[5] *= invSY; te[6] *= invSY;
		te[8] *= invSZ; te[9] *= invSZ; te[10] *= invSZ;

		// Export the rotation.
		quaternion.setFromRotationMatrix(this);

		// Restore the original values.
		te[0] = n00; te[1] = n10; te[2] = n20;
		te[4] = n01; te[5] = n11; te[6] = n21;
		te[8] = n02; te[9] = n12; te[10] = n22;

		// Export the scale.
		scale.x = sx;
		scale.y = sy;
		scale.z = sz;

		return this;

	};

	/**
		 * Creates a perspective matrix.
		 *
		 * @param {Number} left - The distance to the left plane.
		 * @param {Number} right - The distance to the right plane.
		 * @param {Number} top - The distance to the top plane.
		 * @param {Number} bottom - The distance to the bottom plane.
		 * @param {Number} near - The distance to the near plane.
		 * @param {Number} far - The distance to the far plane.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.makePerspective = function makePerspective (left, right, top, bottom, near, far) {

		var te = this.elements;
		var x = 2 * near / (right - left);
		var y = 2 * near / (top - bottom);

		var a = (right + left) / (right - left);
		var b = (top + bottom) / (top - bottom);
		var c = -(far + near) / (far - near);
		var d = -2 * far * near / (far - near);

		te[0] = x; te[4] = 0; te[8] = a; te[12] = 0;
		te[1] = 0; te[5] = y; te[9] = b; te[13] = 0;
		te[2] = 0; te[6] = 0; te[10] = c; te[14] = d;
		te[3] = 0; te[7] = 0; te[11] = -1; te[15] = 0;

		return this;

	};

	/**
		 * Creates an orthographic matrix.
		 *
		 * @param {Number} left - The distance to the left plane.
		 * @param {Number} right - The distance to the right plane.
		 * @param {Number} top - The distance to the top plane.
		 * @param {Number} bottom - The distance to the bottom plane.
		 * @param {Number} near - The distance to the near plane.
		 * @param {Number} far - The distance to the far plane.
		 * @return {Matrix4} This matrix.
		 */

	Matrix4.prototype.makeOrthographic = function makeOrthographic (left, right, top, bottom, near, far) {

		var te = this.elements;
		var w = 1.0 / (right - left);
		var h = 1.0 / (top - bottom);
		var p = 1.0 / (far - near);

		var x = (right + left) * w;
		var y = (top + bottom) * h;
		var z = (far + near) * p;

		te[0] = 2 * w; te[4] = 0; te[8] = 0; te[12] = -x;
		te[1] = 0; te[5] = 2 * h; te[9] = 0; te[13] = -y;
		te[2] = 0; te[6] = 0; te[10] = -2 * p; te[14] = -z;
		te[3] = 0; te[7] = 0; te[11] = 0; te[15] = 1;

		return this;

	};

	/**
		 * Checks if this matrix equals the given one.
		 *
		 * @param {Matrix4} m - A matrix.
		 * @return {Boolean} Whether the matrix are equal.
		 */

	Matrix4.prototype.equals = function equals (matrix) {

		var te = this.elements;
		var me = matrix.elements;

		var result = true;
		var i;

		for(i = 0; result && i < 16; ++i) {

			if(te[i] !== me[i]) {

				result = false;

			}

		}

		return result;

	};

	/**
	 * A vector.
	 *
	 * @type {Vector3}
	 * @private
	 */

	var a$2 = new Vector3();

	/**
	 * A vector.
	 *
	 * @type {Vector3}
	 * @private
	 */

	var b$2 = new Vector3();

	/**
	 * A plane.
	 */

	var Plane = function Plane(normal, constant) {
		if ( normal === void 0 ) normal = new Vector3(1, 0, 0);
		if ( constant === void 0 ) constant = 0;


		/**
			 * The normal.
			 *
			 * @type {Vector3}
			 */

		this.normal = normal;

		/**
			 * The constant.
			 *
			 * @type {Number}
			 */

		this.constant = constant;

	};

	/**
		 * Sets the normal and the constant.
		 *
		 * @param {Vector3} normal - The normal.
		 * @param {Number} constant - The constant.
		 * @return {Plane} This plane.
		 */

	Plane.prototype.set = function set (normal, constant) {

		this.normal.copy(normal);
		this.constant = constant;

		return this;

	};

	/**
		 * Sets the components of this plane.
		 *
		 * @param {Number} x - The X component of the normal.
		 * @param {Number} y - The Y component of the normal.
		 * @param {Number} z - The Z component of the normal.
		 * @param {Number} w - The constant.
		 * @return {Plane} This plane.
		 */

	Plane.prototype.setComponents = function setComponents (x, y, z, w) {

		this.normal.set(x, y, z);
		this.constant = w;

		return this;

	};

	/**
		 * Copies the given plane.
		 *
		 * @param {Plane} p - A plane.
		 * @return {Plane} This plane.
		 */

	Plane.prototype.copy = function copy (p) {

		this.normal.copy(p.normal);
		this.constant = p.constant;

		return this;

	};

	/**
		 * Clones this plane.
		 *
		 * @return {Plane} The cloned plane.
		 */

	Plane.prototype.clone = function clone () {

		return new this.constructor().copy(this);

	};

	/**
		 * Sets this plane from a normal and a coplanar point.
		 *
		 * @param {Vector3} n - The normal.
		 * @param {Vector3} p - The coplanar point.
		 * @return {Sphere} This sphere.
		 */

	Plane.prototype.setFromNormalAndCoplanarPoint = function setFromNormalAndCoplanarPoint (n, p) {

		this.normal.copy(n);
		this.constant = -p.dot(this.normal);

		return this;

	};

	/**
		 * Sets this plane from three distinct coplanar points.
		 *
		 * @param {Vector3} p0 - A coplanar point.
		 * @param {Vector3} p1 - A coplanar point.
		 * @param {Vector3} p2 - A coplanar point.
		 * @return {Plane} This plane.
		 */

	Plane.prototype.setFromCoplanarPoints = function setFromCoplanarPoints (p0, p1, p2) {

		var normal = a$2.subVectors(p2, p1).cross(b$2.subVectors(p0, p1)).normalize();

		this.setFromNormalAndCoplanarPoint(normal, a$2);

		return this;

	};

	/**
		 * Normalizes this plane.
		 *
		 * @return {Plane} This plane.
		 */

	Plane.prototype.normalize = function normalize () {

		var inverseNormalLength = 1.0 / this.normal.length();

		this.normal.multiplyScalar(inverseNormalLength);
		this.constant *= inverseNormalLength;

		return this;

	};

	/**
		 * Negates this plane.
		 *
		 * @return {Plane} This plane.
		 */

	Plane.prototype.negate = function negate () {

		this.normal.negate();
		this.constant = -this.constant;

		return this;

	};

	/**
		 * Calculates the distance from this plane to the given point.
		 *
		 * @param {Vector3} p - A point.
		 * @return {Number} The length.
		 */

	Plane.prototype.distanceToPoint = function distanceToPoint (p) {

		return this.normal.dot(p) + this.constant;

	};

	/**
		 * Calculates the distance from this plane to the given sphere.
		 *
		 * @param {Sphere} s - A sphere.
		 * @return {Number} The length.
		 */

	Plane.prototype.distanceToSphere = function distanceToSphere (s) {

		return this.distanceToPoint(s.center) - s.radius;

	};

	/**
		 * Projects the given point on this plane.
		 *
		 * @param {Vector3} p - A point.
		 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
		 * @return {Vector3} The projected point.
		 */

	Plane.prototype.projectPoint = function projectPoint (p, target) {

		return target.copy(this.normal).multiplyScalar(-this.distanceToPoint(p)).add(p);

	};

	/**
		 * Calculates a coplanar point and returns it.
		 *
		 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
		 * @return {Vector3} A coplanar plane.
		 */

	Plane.prototype.coplanarPoint = function coplanarPoint (target) {

		return target.copy(this.normal).multiplyScalar(-this.constant);

	};

	/**
		 * Translates this plane.
		 *
		 * @param {Vector3} offset - An offset.
		 * @return {Plane} This plane.
		 */

	Plane.prototype.translate = function translate (offset) {

		this.constant -= offset.dot(this.normal);

		return this;

	};

	/**
		 * Finds the point of intersection between this plane and a given line.
		 *
		 * @param {Line3} l - A line.
		 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
		 * @return {Vector3} The intersection point.
		 */

	Plane.prototype.intersectLine = function intersectLine (l, target) {

		var direction = l.delta(a$2);
		var denominator = this.normal.dot(direction);

		if(denominator === 0) {

			// The line is coplanar, return origin.
			if(this.distanceToPoint(l.start) === 0) {

				target.copy(l.start);

			}

		} else {

			var t = -(l.start.dot(this.normal) + this.constant) / denominator;

			if(t >= 0 && t <= 1) {

				target.copy(direction).multiplyScalar(t).add(l.start);

			}

		}

		return target;

	};

	/**
		 * Checks if this plane intersects with the given line.
		 *
		 * @param {Line3} l - A line.
		 * @return {Boolean} Whether this plane intersects with the given line.
		 */

	Plane.prototype.intersectsLine = function intersectsLine (l) {

		var startSign = this.distanceToPoint(l.start);
		var endSign = this.distanceToPoint(l.end);

		return ((startSign < 0 && endSign > 0) || (endSign < 0 && startSign > 0));

	};

	/**
		 * Checks if this plane intersects with the given box.
		 *
		 * @param {Box3} b - A box.
		 * @return {Boolean} Whether this plane intersects with the given box.
		 */

	Plane.prototype.intersectsBox = function intersectsBox (b) {

		return b.intersectsPlane(this);

	};

	/**
		 * Checks if this plane intersects with the given sphere.
		 *
		 * @param {Sphere} s - A sphere.
		 * @return {Boolean} Whether this plane intersects with the given sphere.
		 */

	Plane.prototype.intersectsSphere = function intersectsSphere (s) {

		return s.intersectsPlane(this);

	};

	/**
		 * Checks if this plane equals the given one.
		 *
		 * @param {Plane} v - A plane.
		 * @return {Boolean} Whether this plane equals the given one.
		 */

	Plane.prototype.equals = function equals (p) {

		return (p.normal.equals(this.normal) && (p.constant === this.constant));

	};

	/**
	 * A spherical coordinate system.
	 *
	 * For details see: https://en.wikipedia.org/wiki/Spherical_coordinate_system
	 *
	 * The poles (phi) are at the positive and negative Y-axis. The equator starts
	 * at positive Z.
	 */

	var Spherical = function Spherical(radius, phi, theta) {
		if ( radius === void 0 ) radius = 1;
		if ( phi === void 0 ) phi = 0;
		if ( theta === void 0 ) theta = 0;


		/**
			 * The radius of the sphere.
			 *
			 * @type {Number}
			 */

		this.radius = radius;

		/**
			 * The polar angle, up and down towards the top and bottom pole.
			 *
			 * @type {Number}
			 */

		this.phi = phi;

		/**
			 * The angle around the equator of the sphere.
			 *
			 * @type {Number}
			 */

		this.theta = theta;

	};

	/**
		 * Sets the values of this spherical system.
		 *
		 * @param {Number} radius - The radius.
		 * @param {Number} phi - Phi.
		 * @param {Number} theta - Theta.
		 * @return {Spherical} This spherical system.
		 */

	Spherical.prototype.set = function set (radius, phi, theta) {

		this.radius = radius;
		this.phi = phi;
		this.theta = theta;

		return this;

	};

	/**
		 * Copies the values of the given spherical system.
		 *
		 * @param {Spherical} s - A spherical system.
		 * @return {Spherical} This spherical system.
		 */

	Spherical.prototype.copy = function copy (s) {

		this.radius = s.radius;
		this.phi = s.phi;
		this.theta = s.theta;

		return this;

	};

	/**
		 * Clones this spherical system.
		 *
		 * @return {Spherical} The cloned spherical system.
		 */

	Spherical.prototype.clone = function clone () {

		return new this.constructor().copy(this);

	};

	/**
		 * Sets the values of this spherical system.
		 *
		 * @param {Vector3} v - The vector
		 * @return {Spherical} This spherical system.
		 */

	Spherical.prototype.setFromVector3 = function setFromVector3 (v) {

		this.radius = v.length();

		if(this.radius === 0) {

			this.theta = 0;
			this.phi = 0;

		} else {

			// Calculate the equator angle around the positive Y-axis.
			this.theta = Math.atan2(v.x, v.z);

			// Calculate the polar angle.
			this.phi = Math.acos(Math.min(Math.max(v.y / this.radius, -1), 1));

		}

		return this.makeSafe();

	};

	/**
		 * Restricts phi to [1e-6, PI - 1e-6].
		 *
		 * @return {Spherical} This spherical system.
		 */

	Spherical.prototype.makeSafe = function makeSafe () {

		this.phi = Math.max(1e-6, Math.min(Math.PI - 1e-6, this.phi));

		return this;

	};

	/**
	 * A symmetric 3x3 matrix.
	 */

	var SymmetricMatrix3 = function SymmetricMatrix3() {

		/**
			 * The matrix elements.
			 *
			 * @type {Float32Array}
			 */

		this.elements = new Float32Array([

			1, 0, 0,
			1, 0,
			1

		]);

	};

	/**
		 * Sets the values of this matrix.
		 *
		 * @param {Number} m00 - The value of the first row, first column.
		 * @param {Number} m01 - The value of the first row, second column and the second row, first column.
		 * @param {Number} m02 - The value of the first row, third column and the third row, first column.
		 * @param {Number} m11 - The value of the second row, second column.
		 * @param {Number} m12 - The value of the second row, third column and third row, second column.
		 * @param {Number} m22 - The value of the third row, third column.
		 * @return {SymmetricMatrix3} This matrix.
		 */

	SymmetricMatrix3.prototype.set = function set (m00, m01, m02, m11, m12, m22) {

		var e = this.elements;

		e[0] = m00;
		e[1] = m01; e[3] = m11;
		e[2] = m02; e[4] = m12; e[5] = m22;

		return this;

	};

	/**
		 * Sets this matrix to the identity matrix.
		 *
		 * @return {SymmetricMatrix3} This matrix.
		 */

	SymmetricMatrix3.prototype.identity = function identity () {

		this.set(

			1, 0, 0,
			1, 0,
			1

		);

		return this;

	};

	/**
		 * Copies the values of a given symmetric matrix.
		 *
		 * @param {SymmetricMatrix3} m - A matrix.
		 * @return {SymmetricMatrix3} This matrix.
		 */

	SymmetricMatrix3.prototype.copy = function copy (m) {

		var me = m.elements;

		this.set(

			me[0], me[1], me[2],
			me[3], me[4],
			me[5]

		);

		return this;

	};

	/**
		 * Clones this matrix.
		 *
		 * @return {SymmetricMatrix3} A clone of this matrix.
		 */

	SymmetricMatrix3.prototype.clone = function clone () {

		return new this.constructor().copy(this);

	};

	/**
		 * Copies this symmetric matrix into a given 3x3 matrix.
		 *
		 * @param {Matrix3} m - The target matrix.
		 */

	SymmetricMatrix3.prototype.toMatrix3 = function toMatrix3 (m) {

		var me = m.elements;

		m.set(

			me[0], me[1], me[2],
			me[1], me[3], me[4],
			me[2], me[4], me[5]

		);

	};

	/**
		 * Adds the values of a given symmetric matrix to this one.
		 *
		 * @param {SymmetricMatrix3} m - A matrix.
		 * @return {SymmetricMatrix3} This matrix.
		 */

	SymmetricMatrix3.prototype.add = function add (m) {

		var te = this.elements;
		var me = m.elements;

		te[0] += me[0];
		te[1] += me[1]; te[3] += me[3];
		te[2] += me[2]; te[4] += me[4]; te[5] += me[5];

		return this;

	};

	/**
		 * Calculates the Frobenius norm of this matrix.
		 *
		 * @return {Number} The norm of this matrix.
		 */

	SymmetricMatrix3.prototype.norm = function norm () {

		var e = this.elements;

		var m01m01 = e[1] * e[1];
		var m02m02 = e[2] * e[2];
		var m12m12 = e[4] * e[4];

		return Math.sqrt(

			e[0] * e[0] + m01m01 + m02m02 +
			m01m01 + e[3] * e[3] + m12m12 +
			m02m02 + m12m12 + e[5] * e[5]

		);

	};

	/**
		 * Calculates the absolute sum of all matrix components except for the main
		 * diagonal.
		 *
		 * @return {Number} The offset of this matrix.
		 */

	SymmetricMatrix3.prototype.off = function off () {

		var e = this.elements;

		return Math.sqrt(2 * (

			// Diagonal = [0, 3, 5].
			e[1] * e[1] + e[2] * e[2] + e[4] * e[4]

		));

	};

	/**
		 * Applies this symmetric matrix to a vector.
		 *
		 * @param {Vector3} v - The vector to modify.
		 * @return {Vector3} The modified vector.
		 */

	SymmetricMatrix3.prototype.applyToVector3 = function applyToVector3 (v) {

		var x = v.x, y = v.y, z = v.z;
		var e = this.elements;

		v.x = e[0] * x + e[1] * y + e[2] * z;
		v.y = e[1] * x + e[3] * y + e[4] * z;
		v.z = e[2] * x + e[4] * y + e[5] * z;

		return v;

	};

	/**
		 * Checks if this matrix equals the given one.
		 *
		 * @param {SymmetricMatrix3} m - A matrix.
		 * @return {Boolean} Whether the matrix are equal.
		 */

	SymmetricMatrix3.prototype.equals = function equals (matrix) {

		var te = this.elements;
		var me = matrix.elements;

		var result = true;
		var i;

		for(i = 0; result && i < 6; ++i) {

			if(te[i] !== me[i]) {

				result = false;

			}

		}

		return result;

	};

	/**
		 * Calculates the linear index of an element from this matrix.
		 *
		 * Let N be the dimension of the symmetric matrix:
		 *
		 *     index = N * (N - 1) / 2 - (N - i) * (N - i - 1) / 2 + j
		 *
		 * @param {Number} i - The row.
		 * @param {Number} j - The column.
		 * @return {Number} The index into the elements of this matrix.
		 */

	SymmetricMatrix3.calculateIndex = function calculateIndex (i, j) {

		return (3 - (3 - i) * (2 - i) / 2 + j);

	};

	/**
	 * A vector with four components.
	 */

	var Vector4 = function Vector4(x, y, z, w) {
		if ( x === void 0 ) x = 0;
		if ( y === void 0 ) y = 0;
		if ( z === void 0 ) z = 0;
		if ( w === void 0 ) w = 0;


		/**
			 * The X component.
			 *
			 * @type {Number}
			 */

		this.x = x;

		/**
			 * The Y component.
			 *
			 * @type {Number}
			 */

		this.y = y;

		/**
			 * The Z component.
			 *
			 * @type {Number}
			 */

		this.z = z;

		/**
			 * The W component.
			 *
			 * @type {Number}
			 */

		this.w = w;

	};

	/**
		 * Sets the values of this vector
		 *
		 * @param {Number} x - The X component.
		 * @param {Number} y - The Y component.
		 * @param {Number} z - The Z component.
		 * @param {Number} w - The W component.
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.set = function set (x, y, z, w) {

		this.x = x;
		this.y = y;
		this.z = z;
		this.w = w;

		return this;

	};

	/**
		 * Copies the values of another vector.
		 *
		 * @param {Vector4} v - A vector.
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.copy = function copy (v) {

		this.x = v.x;
		this.y = v.y;
		this.z = v.z;
		this.w = v.w;

		return this;

	};

	/**
		 * Clones this vector.
		 *
		 * @return {Vector4} A clone of this vector.
		 */

	Vector4.prototype.clone = function clone () {

		return new this.constructor(this.x, this.y, this.z, this.w);

	};

	/**
		 * Copies values from an array.
		 *
		 * @param {Number[]} array - An array.
		 * @param {Number} offset - An offset.
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.fromArray = function fromArray (array, offset) {
			if ( offset === void 0 ) offset = 0;


		this.x = array[offset];
		this.y = array[offset + 1];
		this.z = array[offset + 2];
		this.w = array[offset + 3];

		return this;

	};

	/**
		 * Stores this vector in an array.
		 *
		 * @param {Array} [array] - A target array.
		 * @param {Number} offset - An offset.
		 * @return {Number[]} The array.
		 */

	Vector4.prototype.toArray = function toArray (array, offset) {
			if ( array === void 0 ) array = [];
			if ( offset === void 0 ) offset = 0;


		array[offset] = this.x;
		array[offset + 1] = this.y;
		array[offset + 2] = this.z;
		array[offset + 3] = this.w;

		return array;

	};

	/**
		 * Stores the axis angle from the given quaternion in this vector.
		 *
		 * For more details see:
		 *  http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToAngle/index.htm
		 *
		 * @param {Quaternion} q - A quaternion. Assumed to be normalized
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.setAxisAngleFromQuaternion = function setAxisAngleFromQuaternion (q) {

		this.w = 2 * Math.acos(q.w);

		var s = Math.sqrt(1 - q.w * q.w);

		if(s < 1e-4) {

			this.x = 1;
			this.y = 0;
			this.z = 0;

		} else {

			this.x = q.x / s;
			this.y = q.y / s;
			this.z = q.z / s;

		}

		return this;

	};

	/**
		 * Stores the axis angle from the given rotation matrix in this vector.
		 *
		 * For more details see:
		 *  http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToAngle/index.htm
		 *
		 * @param {Matrix4} m - A matrix. The upper 3x3 must be a pure rotation matrix (i.e. unscaled).
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.setAxisAngleFromRotationMatrix = function setAxisAngleFromRotationMatrix (m) {

		// Margin to allow for rounding errors.
		var E = 0.01;
		// Margin to distinguish between 0 and 180 degrees.
		var H = 0.1;

		var me = m.elements;
		var m00 = me[0], m01 = me[4], m02 = me[8];
		var m10 = me[1], m11 = me[5], m12 = me[9];
		var m20 = me[2], m21 = me[6], m22 = me[10];

		var angle;
		var x, y, z;
		var xx, yy, zz;
		var xy, xz, yz;
		var s;

		if((Math.abs(m01 - m10) < E) && (Math.abs(m02 - m20) < E) && (Math.abs(m12 - m21) < E)) {

			/* Singularity found. First, check for identity matrix which must have +1
			for all terms in the leading diagonal and zero in other terms. */
			if((Math.abs(m01 + m10) < H) && (Math.abs(m02 + m20) < H) && (Math.abs(m12 + m21) < H) && (Math.abs(m00 + m11 + m22 - 3) < H)) {

				// This singularity is the identity matrix. The angle is zero.
				this.set(1, 0, 0, 0);

			} else {

				// The angle is 180.
				angle = Math.PI;

				xx = (m00 + 1) / 2;
				yy = (m11 + 1) / 2;
				zz = (m22 + 1) / 2;
				xy = (m01 + m10) / 4;
				xz = (m02 + m20) / 4;
				yz = (m12 + m21) / 4;

				if((xx > yy) && (xx > zz)) {

					// m00 is the largest diagonal term.
					if(xx < E) {

						x = 0;
						y = 0.707106781;
						z = 0.707106781;

					} else {

						x = Math.sqrt(xx);
						y = xy / x;
						z = xz / x;

					}

				} else if(yy > zz) {

					// m11 is the largest diagonal term.
					if(yy < E) {

						x = 0.707106781;
						y = 0;
						z = 0.707106781;

					} else {

						y = Math.sqrt(yy);
						x = xy / y;
						z = yz / y;

					}

				} else {

					// m22 is the largest diagonal term.
					if(zz < E) {

						x = 0.707106781;
						y = 0.707106781;
						z = 0;

					} else {

						z = Math.sqrt(zz);
						x = xz / z;
						y = yz / z;

					}

				}

				this.set(x, y, z, angle);

			}

		} else {

			// There are no singularities.
			s = Math.sqrt(
				(m21 - m12) * (m21 - m12) +
				(m02 - m20) * (m02 - m20) +
				(m10 - m01) * (m10 - m01)
			);

			// Prevent division by zero.
			if(Math.abs(s) < 0.001) { s = 1; }

			this.x = (m21 - m12) / s;
			this.y = (m02 - m20) / s;
			this.z = (m10 - m01) / s;
			this.w = Math.acos((m00 + m11 + m22 - 1) / 2);

		}

		return this;

	};

	/**
		 * Adds a vector to this one.
		 *
		 * @param {Vector4} v - The vector to add.
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.add = function add (v) {

		this.x += v.x;
		this.y += v.y;
		this.z += v.z;
		this.w += v.w;

		return this;

	};

	/**
		 * Adds a scalar to this vector.
		 *
		 * @param {Number} s - The scalar to add.
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.addScalar = function addScalar (s) {

		this.x += s;
		this.y += s;
		this.z += s;
		this.w += s;

		return this;

	};

	/**
		 * Sets this vector to the sum of two given vectors.
		 *
		 * @param {Vector4} a - A vector.
		 * @param {Vector4} b - Another vector.
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.addVectors = function addVectors (a, b) {

		this.x = a.x + b.x;
		this.y = a.y + b.y;
		this.z = a.z + b.z;
		this.w = a.w + b.w;

		return this;

	};

	/**
		 * Adds a scaled vector to this one.
		 *
		 * @param {Vector4} v - The vector to scale and add.
		 * @param {Number} s - A scalar.
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.addScaledVector = function addScaledVector (v, s) {

		this.x += v.x * s;
		this.y += v.y * s;
		this.z += v.z * s;
		this.w += v.w * s;

		return this;

	};

	/**
		 * Subtracts a vector from this vector.
		 *
		 * @param {Vector4} v - The vector to subtract.
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.sub = function sub (v) {

		this.x -= v.x;
		this.y -= v.y;
		this.z -= v.z;
		this.w -= v.w;

		return this;

	};

	/**
		 * Subtracts a scalar to this vector.
		 *
		 * @param {Number} s - The scalar to subtract.
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.subScalar = function subScalar (s) {

		this.x -= s;
		this.y -= s;
		this.z -= s;
		this.w -= s;

		return this;

	};

	/**
		 * Sets this vector to the difference between two given vectors.
		 *
		 * @param {Vector4} a - A vector.
		 * @param {Vector4} b - A second vector.
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.subVectors = function subVectors (a, b) {

		this.x = a.x - b.x;
		this.y = a.y - b.y;
		this.z = a.z - b.z;
		this.w = a.w - b.w;

		return this;

	};

	/**
		 * Multiplies this vector with another vector.
		 *
		 * @param {Vector4} v - A vector.
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.multiply = function multiply (v) {

		this.x *= v.x;
		this.y *= v.y;
		this.z *= v.z;
		this.w *= v.w;

		return this;

	};

	/**
		 * Multiplies this vector with a given scalar.
		 *
		 * @param {Number} s - A scalar.
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.multiplyScalar = function multiplyScalar (s) {

		this.x *= s;
		this.y *= s;
		this.z *= s;
		this.w *= s;

		return this;

	};

	/**
		 * Sets this vector to the product of two given vectors.
		 *
		 * @param {Vector4} a - A vector.
		 * @param {Vector4} b - Another vector.
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.multiplyVectors = function multiplyVectors (a, b) {

		this.x = a.x * b.x;
		this.y = a.y * b.y;
		this.z = a.z * b.z;
		this.w = a.w * b.w;

		return this;

	};

	/**
		 * Divides this vector by another vector.
		 *
		 * @param {Vector4} v - A vector.
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.divide = function divide (v) {

		this.x /= v.x;
		this.y /= v.y;
		this.z /= v.z;
		this.w /= v.w;

		return this;

	};

	/**
		 * Divides this vector by a given scalar.
		 *
		 * @param {Number} s - A scalar.
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.divideScalar = function divideScalar (s) {

		this.x /= s;
		this.y /= s;
		this.z /= s;
		this.w /= s;

		return this;

	};

	/**
		 * Applies a matrix to this vector.
		 *
		 * @param {Matrix4} m - A matrix.
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.applyMatrix4 = function applyMatrix4 (m) {

		var x = this.x, y = this.y, z = this.z, w = this.w;
		var e = m.elements;

		this.x = e[0] * x + e[4] * y + e[8] * z + e[12] * w;
		this.y = e[1] * x + e[5] * y + e[9] * z + e[13] * w;
		this.z = e[2] * x + e[6] * y + e[10] * z + e[14] * w;
		this.w = e[3] * x + e[7] * y + e[11] * z + e[15] * w;

		return this;

	};

	/**
		 * Negates this vector.
		 *
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.negate = function negate () {

		this.x = -this.x;
		this.y = -this.y;
		this.z = -this.z;
		this.w = -this.w;

		return this;

	};

	/**
		 * Calculates the dot product with another vector.
		 *
		 * @param {Vector4} v - A vector.
		 * @return {Number} The dot product.
		 */

	Vector4.prototype.dot = function dot (v) {

		return this.x * v.x + this.y * v.y + this.z * v.z + this.w * v.w;

	};

	/**
		 * Calculates the Manhattan length of this vector.
		 *
		 * @return {Number} The length.
		 */

	Vector4.prototype.lengthManhattan = function lengthManhattan () {

		return Math.abs(this.x) + Math.abs(this.y) + Math.abs(this.z) + Math.abs(this.w);

	};

	/**
		 * Calculates the squared length of this vector.
		 *
		 * @return {Number} The squared length.
		 */

	Vector4.prototype.lengthSquared = function lengthSquared () {

		return (
			this.x * this.x +
			this.y * this.y +
			this.z * this.z +
			this.w * this.w
		);

	};

	/**
		 * Calculates the length of this vector.
		 *
		 * @return {Number} The length.
		 */

	Vector4.prototype.length = function length () {

		return Math.sqrt(
			this.x * this.x +
			this.y * this.y +
			this.z * this.z +
			this.w * this.w
		);

	};

	/**
		 * Calculates the Manhattan distance to a given vector.
		 *
		 * @param {Vector4} v - A vector.
		 * @return {Number} The distance.
		 */

	Vector4.prototype.distanceToManhattan = function distanceToManhattan (v) {

		return (
			Math.abs(this.x - v.x) +
			Math.abs(this.y - v.y) +
			Math.abs(this.z - v.z) +
			Math.abs(this.w - v.w)
		);

	};

	/**
		 * Calculates the squared distance to a given vector.
		 *
		 * @param {Vector4} v - A vector.
		 * @return {Number} The squared distance.
		 */

	Vector4.prototype.distanceToSquared = function distanceToSquared (v) {

		var dx = this.x - v.x;
		var dy = this.y - v.y;
		var dz = this.z - v.z;
		var dw = this.w - v.w;

		return dx * dx + dy * dy + dz * dz + dw * dw;

	};

	/**
		 * Calculates the distance to a given vector.
		 *
		 * @param {Vector4} v - A vector.
		 * @return {Number} The distance.
		 */

	Vector4.prototype.distanceTo = function distanceTo (v) {

		return Math.sqrt(this.distanceToSquared(v));

	};

	/**
		 * Normalizes this vector.
		 *
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.normalize = function normalize () {

		return this.divideScalar(this.length());

	};

	/**
		 * Sets the length of this vector.
		 *
		 * @param {Number} length - The new length.
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.setLength = function setLength (length) {

		return this.normalize().multiplyScalar(length);

	};

	/**
		 * Adopts the min value for each component of this vector and the given one.
		 *
		 * @param {Vector4} v - A vector.
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.min = function min (v) {

		this.x = Math.min(this.x, v.x);
		this.y = Math.min(this.y, v.y);
		this.z = Math.min(this.z, v.z);
		this.w = Math.min(this.w, v.w);

		return this;

	};

	/**
		 * Adopts the max value for each component of this vector and the given one.
		 *
		 * @param {Vector4} v - A vector.
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.max = function max (v) {

		this.x = Math.max(this.x, v.x);
		this.y = Math.max(this.y, v.y);
		this.z = Math.max(this.z, v.z);
		this.w = Math.max(this.w, v.w);

		return this;

	};

	/**
		 * Clamps this vector.
		 *
		 * @param {Vector4} min - The lower bounds. Assumed to be smaller than max.
		 * @param {Vector4} max - The upper bounds. Assumed to be greater than min.
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.clamp = function clamp (min, max) {

		this.x = Math.max(min.x, Math.min(max.x, this.x));
		this.y = Math.max(min.y, Math.min(max.y, this.y));
		this.z = Math.max(min.z, Math.min(max.z, this.z));
		this.w = Math.max(min.w, Math.min(max.w, this.w));

		return this;

	};

	/**
		 * Floors this vector.
		 *
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.floor = function floor () {

		this.x = Math.floor(this.x);
		this.y = Math.floor(this.y);
		this.z = Math.floor(this.z);
		this.w = Math.floor(this.w);

		return this;

	};

	/**
		 * Ceils this vector.
		 *
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.ceil = function ceil () {

		this.x = Math.ceil(this.x);
		this.y = Math.ceil(this.y);
		this.z = Math.ceil(this.z);
		this.w = Math.ceil(this.w);

		return this;

	};

	/**
		 * Rounds this vector.
		 *
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.round = function round () {

		this.x = Math.round(this.x);
		this.y = Math.round(this.y);
		this.z = Math.round(this.z);
		this.w = Math.round(this.w);

		return this;

	};

	/**
		 * Lerps towards the given vector.
		 *
		 * @param {Vector4} v - The target vector.
		 * @param {Number} alpha - The lerp factor.
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.lerp = function lerp (v, alpha) {

		this.x += (v.x - this.x) * alpha;
		this.y += (v.y - this.y) * alpha;
		this.z += (v.z - this.z) * alpha;
		this.w += (v.w - this.w) * alpha;

		return this;

	};

	/**
		 * Sets this vector to the lerp result of the given vectors.
		 *
		 * @param {Vector4} v1 - A base vector.
		 * @param {Vector4} v2 - The target vector.
		 * @param {Number} alpha - The lerp factor.
		 * @return {Vector4} This vector.
		 */

	Vector4.prototype.lerpVectors = function lerpVectors (v1, v2, alpha) {

		return this.subVectors(v2, v1).multiplyScalar(alpha).add(v1);

	};

	/**
		 * Checks if this vector equals the given one.
		 *
		 * @param {Vector4} v - A vector.
		 * @return {Boolean} Whether this vector equals the given one.
		 */

	Vector4.prototype.equals = function equals (v) {

		return (v.x === this.x && v.y === this.y && v.z === this.z && v.w === this.w);

	};

	/**
	 * Mathematical data structures.
	 *
	 * @module math-ds
	 */

	/**
	 * A vector.
	 *
	 * @type {Vector3}
	 * @private
	 */

	var c$2 = new Vector3();

	/**
	 * An octant.
	 */

	var Octant = function Octant(min, max) {
		if ( min === void 0 ) min = new Vector3();
		if ( max === void 0 ) max = new Vector3();


		/**
			 * The lower bounds of this octant.
			 *
			 * @type {Vector3}
			 */

		this.min = min;

		/**
			 * The upper bounds of the octant.
			 *
			 * @type {Vector3}
			 */

		this.max = max;

		/**
			 * The children of this octant.
			 *
			 * @type {Octant[]}
			 * @default null
			 */

		this.children = null;

	};

	/**
		 * Computes the center of this octant.
		 *
		 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
		 * @return {Vector3} A vector that describes the center of this octant.
		 */

	Octant.prototype.getCenter = function getCenter (target) {
			if ( target === void 0 ) target = new Vector3();


		return target.addVectors(this.min, this.max).multiplyScalar(0.5);

	};

	/**
		 * Computes the size of this octant.
		 *
		 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
		 * @return {Vector3} A vector that describes the size of this octant.
		 */

	Octant.prototype.getDimensions = function getDimensions (target) {
			if ( target === void 0 ) target = new Vector3();


		return target.subVectors(this.max, this.min);

	};

	/**
		 * Splits this octant into eight smaller ones.
		 */

	Octant.prototype.split = function split () {
			var this$1 = this;


		var min = this.min;
		var max = this.max;
		var mid = this.getCenter(c$2);

		var children = this.children = [

			null, null,
			null, null,
			null, null,
			null, null

		];

		var i, combination;

		for(i = 0; i < 8; ++i) {

			combination = pattern[i];

			children[i] = new this$1.constructor(

				new Vector3(
					(combination[0] === 0) ? min.x : mid.x,
					(combination[1] === 0) ? min.y : mid.y,
					(combination[2] === 0) ? min.z : mid.z
				),

				new Vector3(
					(combination[0] === 0) ? mid.x : max.x,
					(combination[1] === 0) ? mid.y : max.y,
					(combination[2] === 0) ? mid.z : max.z
				)

			);

		}

	};

	/**
	 * A binary pattern that describes the standard octant layout:
	 *
	 * ```text
	 *    3____7
	 *  2/___6/|
	 *  | 1__|_5
	 *  0/___4/
	 * ```
	 *
	 * This common layout is crucial for positional assumptions.
	 *
	 * @type {Uint8Array[]}
	 */

	var pattern = [

		new Uint8Array([0, 0, 0]),
		new Uint8Array([0, 0, 1]),
		new Uint8Array([0, 1, 0]),
		new Uint8Array([0, 1, 1]),

		new Uint8Array([1, 0, 0]),
		new Uint8Array([1, 0, 1]),
		new Uint8Array([1, 1, 0]),
		new Uint8Array([1, 1, 1])

	];

	/**
	 * Describes all possible octant corner connections.
	 *
	 * @type {Uint8Array[]}
	 */

	var edges = [

		// X-Axis.
		new Uint8Array([0, 4]),
		new Uint8Array([1, 5]),
		new Uint8Array([2, 6]),
		new Uint8Array([3, 7]),

		// Y-Axis.
		new Uint8Array([0, 2]),
		new Uint8Array([1, 3]),
		new Uint8Array([4, 6]),
		new Uint8Array([5, 7]),

		// Z-Axis.
		new Uint8Array([0, 1]),
		new Uint8Array([2, 3]),
		new Uint8Array([4, 5]),
		new Uint8Array([6, 7])

	];

	/**
	 * A vector.
	 *
	 * @type {Vector3}
	 * @private
	 */

	var c = new Vector3();

	/**
	 * A cubic octant.
	 */

	var CubicOctant = function CubicOctant(min, size) {
		if ( min === void 0 ) min = new Vector3();
		if ( size === void 0 ) size = 0;


		/**
			 * The lower bounds of this octant.
			 *
			 * @type {Vector3}
			 */

		this.min = min;

		/**
			 * The size of this octant.
			 *
			 * @type {Number}
			 */

		this.size = size;

		/**
			 * The children of this octant.
			 *
			 * @type {CubicOctant[]}
			 * @default null
			 */

		this.children = null;

	};

	var prototypeAccessors = { max: { configurable: true } };

	/**
		 * The upper bounds of this octant.
		 *
		 * Accessing this property always creates a new vector.
		 *
		 * @type {Vector3}
		 */

	prototypeAccessors.max.get = function () { return this.min.clone().addScalar(this.size); };

	/**
		 * Computes the center of this octant.
		 *
		 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
		 * @return {Vector3} A vector that describes the center of this octant.
		 */

	CubicOctant.prototype.getCenter = function getCenter (target) {
			if ( target === void 0 ) target = new Vector3();


		return target.copy(this.min).addScalar(this.size * 0.5);

	};

	/**
		 * Returns the size of this octant as a vector.
		 *
		 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
		 * @return {Vector3} A vector that describes the size of this octant.
		 */

	CubicOctant.prototype.getDimensions = function getDimensions (target) {
			if ( target === void 0 ) target = new Vector3();


		return target.set(this.size, this.size, this.size);

	};

	/**
		 * Splits this octant into eight smaller ones.
		 */

	CubicOctant.prototype.split = function split () {
			var this$1 = this;


		var min = this.min;
		var mid = this.getCenter(c);
		var halfSize = this.size * 0.5;

		var children = this.children = [

			null, null,
			null, null,
			null, null,
			null, null

		];

		var i, combination;

		for(i = 0; i < 8; ++i) {

			combination = pattern[i];

			children[i] = new this$1.constructor(

				new Vector3(
					(combination[0] === 0) ? min.x : mid.x,
					(combination[1] === 0) ? min.y : mid.y,
					(combination[2] === 0) ? min.z : mid.z
				),

				halfSize

			);

		}

	};

	Object.defineProperties( CubicOctant.prototype, prototypeAccessors );

	/**
	 * A basic iterator result.
	 *
	 * The next method of an iterator always has to return an object with
	 * appropriate properties including done and value.
	 */

	var IteratorResult = function IteratorResult(value, done) {
		if ( value === void 0 ) value = null;
		if ( done === void 0 ) done = false;


		/**
			 * An arbitrary value returned by the iterator.
			 *
			 * @type Object
			 * @default null
			 */

		this.value = value;

		/**
			 * Whether this result is past the end of the iterated sequence.
			 *
			 * @type Boolean
			 * @default false
			 */

		this.done = done;

	};

	/**
		 * Resets this iterator result.
		 */

	IteratorResult.prototype.reset = function reset () {

		this.value = null;
		this.done = false;

	};

	/**
	 * A compilation of the library components.
	 *
	 * @module iterator-result
	 */

	/**
	 * A 3D box.
	 *
	 * @type {Box3}
	 * @private
	 */

	var b$4 = new Box3();

	/**
	 * An octant iterator.
	 *
	 * @implements {Iterator}
	 * @implements {Iterable}
	 */

	var OctantIterator = function OctantIterator(octree, region) {
		if ( region === void 0 ) region = null;


		/**
			 * The octree.
			 *
			 * @type {Octree}
			 * @private
			 */

		this.octree = octree;

		/**
			 * A region used for octree culling.
			 *
			 * @type {Frustum|Box3}
			 * @default null
			 */

		this.region = region;

		/**
			 * Whether this iterator should respect the cull region.
			 *
			 * @type {Boolean}
			 * @default false
			 */

		this.cull = (region !== null);

		/**
			 * An iterator result.
			 *
			 * @type {IteratorResult}
			 * @private
			 */

		this.result = new IteratorResult();

		/**
			 * An octant trace.
			 *
			 * @type {Octant[]}
			 * @private
			 */

		this.trace = null;

		/**
			 * Iteration indices.
			 *
			 * @type {Number[]}
			 * @private
			 */

		this.indices = null;

		this.reset();

	};

	/**
		 * Resets this iterator.
		 *
		 * @return {OctantIterator} This iterator.
		 */

	OctantIterator.prototype.reset = function reset () {

		var root = this.octree.root;

		this.trace = [];
		this.indices = [];

		if(root !== null) {

			b$4.min = root.min;
			b$4.max = root.max;

			if(!this.cull || this.region.intersectsBox(b$4)) {

				this.trace.push(root);
				this.indices.push(0);

			}

		}

		this.result.reset();

		return this;

	};

	/**
		 * Iterates over the leaf octants.
		 *
		 * @return {IteratorResult} The next leaf octant.
		 */

	OctantIterator.prototype.next = function next () {

		var cull = this.cull;
		var region = this.region;
		var indices = this.indices;
		var trace = this.trace;

		var octant = null;
		var depth = trace.length - 1;

		var index, children, child;

		while(octant === null && depth >= 0) {

			index = indices[depth];
			children = trace[depth].children;

			++indices[depth];

			if(index < 8) {

				if(children !== null) {

					child = children[index];

					if(cull) {

						b$4.min = child.min;
						b$4.max = child.max;

						if(!region.intersectsBox(b$4)) {

							// Cull this octant.
							continue;

						}

					}

					trace.push(child);
					indices.push(0);

					++depth;

				} else {

					octant = trace.pop();
					indices.pop();

				}

			} else {

				trace.pop();
				indices.pop();

				--depth;

			}

		}

		this.result.value = octant;
		this.result.done = (octant === null);

		return this.result;

	};

	/**
		 * Called when this iterator will no longer be run to completion.
		 *
		 * @param {Object} value - An interator result value.
		 * @return {IteratorResult} - A premature completion result.
		 */

	OctantIterator.prototype.return = function return$1 (value) {

		this.result.value = value;
		this.result.done = true;

		return this.result;

	};

	/**
		 * Returns this iterator.
		 *
		 * @return {OctantIterator} An iterator.
		 */

	OctantIterator.prototype[Symbol.iterator] = function () {

		return this;

	};

	/**
	 * Contains bytes used for bitwise operations. The last byte is used to store
	 * raycasting flags.
	 *
	 * @type Uint8Array
	 * @private
	 */

	var flags = new Uint8Array([0, 1, 2, 3, 4, 5, 6, 7, 0]);

	/**
	 * A lookup-table containing octant ids. Used to determine the exit plane from
	 * an octant.
	 *
	 * @type {Uint8Array[]}
	 * @private
	 */

	var octantTable = [

		new Uint8Array([4, 2, 1]),
		new Uint8Array([5, 3, 8]),
		new Uint8Array([6, 8, 3]),
		new Uint8Array([7, 8, 8]),
		new Uint8Array([8, 6, 5]),
		new Uint8Array([8, 7, 8]),
		new Uint8Array([8, 8, 7]),
		new Uint8Array([8, 8, 8])

	];

	/**
	 * Finds the entry plane of the first octant that a ray travels through.
	 *
	 * Determining the first octant requires knowing which of the t0s is the
	 * largest. The tms of the other axes must also be compared against that
	 * largest t0.
	 *
	 * @private
	 * @param {Number} tx0 - Ray projection parameter.
	 * @param {Number} ty0 - Ray projection parameter.
	 * @param {Number} tz0 - Ray projection parameter.
	 * @param {Number} txm - Ray projection parameter mean.
	 * @param {Number} tym - Ray projection parameter mean.
	 * @param {Number} tzm - Ray projection parameter mean.
	 * @return {Number} The index of the first octant that the ray travels through.
	 */

	function findEntryOctant(tx0, ty0, tz0, txm, tym, tzm) {

		var entry = 0;

		// Find the entry plane.
		if(tx0 > ty0 && tx0 > tz0) {

			// YZ-plane.
			if(tym < tx0) { entry |= 2; }
			if(tzm < tx0) { entry |= 1; }

		} else if(ty0 > tz0) {

			// XZ-plane.
			if(txm < ty0) { entry |= 4; }
			if(tzm < ty0) { entry |= 1; }

		} else {

			// XY-plane.
			if(txm < tz0) { entry |= 4; }
			if(tym < tz0) { entry |= 2; }

		}

		return entry;

	}

	/**
	 * Finds the next octant that intersects with the ray based on the exit plane of
	 * the current one.
	 *
	 * @private
	 * @param {Number} currentOctant - The index of the current octant.
	 * @param {Number} tx1 - Ray projection parameter.
	 * @param {Number} ty1 - Ray projection parameter.
	 * @param {Number} tz1 - Ray projection parameter.
	 * @return {Number} The index of the next octant that the ray travels through.
	 */

	function findNextOctant(currentOctant, tx1, ty1, tz1) {

		var min;
		var exit = 0;

		// Find the exit plane.
		if(tx1 < ty1) {

			min = tx1;
			exit = 0; // YZ-plane.

		} else {

			min = ty1;
			exit = 1; // XZ-plane.

		}

		if(tz1 < min) {

			exit = 2; // XY-plane.

		}

		return octantTable[currentOctant][exit];

	}

	/**
	 * Finds all octants that intersect with the given ray.
	 *
	 * @private
	 * @param {Octant} octant - The current octant.
	 * @param {Number} tx0 - Ray projection parameter. Initial tx0 = (minX - rayOriginX) / rayDirectionX.
	 * @param {Number} ty0 - Ray projection parameter. Initial ty0 = (minY - rayOriginY) / rayDirectionY.
	 * @param {Number} tz0 - Ray projection parameter. Initial tz0 = (minZ - rayOriginZ) / rayDirectionZ.
	 * @param {Number} tx1 - Ray projection parameter. Initial tx1 = (maxX - rayOriginX) / rayDirectionX.
	 * @param {Number} ty1 - Ray projection parameter. Initial ty1 = (maxY - rayOriginY) / rayDirectionY.
	 * @param {Number} tz1 - Ray projection parameter. Initial tz1 = (maxZ - rayOriginZ) / rayDirectionZ.
	 * @param {Raycaster} raycaster - The raycaster.
	 * @param {Array} intersects - An array to be filled with the intersecting octants.
	 */

	function raycastOctant(octant, tx0, ty0, tz0, tx1, ty1, tz1, raycaster, intersects) {

		var children = octant.children;

		var currentOctant;
		var txm, tym, tzm;

		if(tx1 >= 0.0 && ty1 >= 0.0 && tz1 >= 0.0) {

			if(children === null) {

				// Leaf.
				intersects.push(octant);

			} else {

				// Compute means.
				txm = 0.5 * (tx0 + tx1);
				tym = 0.5 * (ty0 + ty1);
				tzm = 0.5 * (tz0 + tz1);

				currentOctant = findEntryOctant(tx0, ty0, tz0, txm, tym, tzm);

				do {

					/* The possibilities for the next node are passed in the same respective
					 * order as the t-values. Hence, if the first value is found to be the
					 * greatest, the fourth one will be returned. If the second value is the
					 * greatest, the fifth one will be returned, etc.
					 */

					switch(currentOctant) {

						case 0:
							raycastOctant(children[flags[8]], tx0, ty0, tz0, txm, tym, tzm, raycaster, intersects);
							currentOctant = findNextOctant(currentOctant, txm, tym, tzm);
							break;

						case 1:
							raycastOctant(children[flags[8] ^ flags[1]], tx0, ty0, tzm, txm, tym, tz1, raycaster, intersects);
							currentOctant = findNextOctant(currentOctant, txm, tym, tz1);
							break;

						case 2:
							raycastOctant(children[flags[8] ^ flags[2]], tx0, tym, tz0, txm, ty1, tzm, raycaster, intersects);
							currentOctant = findNextOctant(currentOctant, txm, ty1, tzm);
							break;

						case 3:
							raycastOctant(children[flags[8] ^ flags[3]], tx0, tym, tzm, txm, ty1, tz1, raycaster, intersects);
							currentOctant = findNextOctant(currentOctant, txm, ty1, tz1);
							break;

						case 4:
							raycastOctant(children[flags[8] ^ flags[4]], txm, ty0, tz0, tx1, tym, tzm, raycaster, intersects);
							currentOctant = findNextOctant(currentOctant, tx1, tym, tzm);
							break;

						case 5:
							raycastOctant(children[flags[8] ^ flags[5]], txm, ty0, tzm, tx1, tym, tz1, raycaster, intersects);
							currentOctant = findNextOctant(currentOctant, tx1, tym, tz1);
							break;

						case 6:
							raycastOctant(children[flags[8] ^ flags[6]], txm, tym, tz0, tx1, ty1, tzm, raycaster, intersects);
							currentOctant = findNextOctant(currentOctant, tx1, ty1, tzm);
							break;

						case 7:
							raycastOctant(children[flags[8] ^ flags[7]], txm, tym, tzm, tx1, ty1, tz1, raycaster, intersects);
							// Far top right octant. No other octants can be reached from here.
							currentOctant = 8;
							break;

					}

				} while(currentOctant < 8);

			}

		}

	}

	/**
	 * The dimensions of an octree.
	 *
	 * @type {Vector3}
	 * @private
	 */

	var dimensions = new Vector3();

	/**
	 * The half dimensions of an octree.
	 *
	 * @type {Vector3}
	 * @private
	 */

	var halfDimensions = new Vector3();

	/**
	 * The center of an octree.
	 *
	 * @type {Vector3}
	 * @private
	 */

	var center = new Vector3();

	/**
	 * The lower bounds of an octree.
	 *
	 * @type {Vector3}
	 * @private
	 */

	var min = new Vector3();

	/**
	 * The upper bounds of an octree.
	 *
	 * @type {Vector3}
	 * @private
	 */

	var max = new Vector3();

	/**
	 * A ray direction.
	 *
	 * @type {Vector3}
	 * @private
	 */

	var direction = new Vector3();

	/**
	 * A ray origin.
	 *
	 * @type {Vector3}
	 * @private
	 */

	var origin = new Vector3();

	/**
	 * An octree raycaster.
	 *
	 * Based on:
	 *  "An Efficient Parametric Algorithm for Octree Traversal"
	 *  by J. Revelles et al. (2000).
	 */

	var OctreeRaycaster = function OctreeRaycaster () {};

	OctreeRaycaster.intersectOctree = function intersectOctree (octree, raycaster, intersects) {

		var invDirX, invDirY, invDirZ;
		var tx0, tx1, ty0, ty1, tz0, tz1;

		octree.getDimensions(dimensions);
		halfDimensions.copy(dimensions).multiplyScalar(0.5);

		// Translate the octree extents to the center of the octree.
		min.copy(octree.min).sub(octree.min);
		max.copy(octree.max).sub(octree.min);

		direction.copy(raycaster.ray.direction);
		origin.copy(raycaster.ray.origin);

		// Translate the ray to the center of the octree.
		origin.sub(octree.getCenter(center)).add(halfDimensions);

		// Reset all flags.
		flags[8] = flags[0];

		// Handle rays with negative directions.
		if(direction.x < 0.0) {

			origin.x = dimensions.x - origin.x;
			direction.x = -direction.x;
			flags[8] |= flags[4];

		}

		if(direction.y < 0.0) {

			origin.y = dimensions.y - origin.y;
			direction.y = -direction.y;
			flags[8] |= flags[2];

		}

		if(direction.z < 0.0) {

			origin.z = dimensions.z - origin.z;
			direction.z = -direction.z;
			flags[8] |= flags[1];

		}

		// Improve IEEE double stability.
		invDirX = 1.0 / direction.x;
		invDirY = 1.0 / direction.y;
		invDirZ = 1.0 / direction.z;

		// Project the ray to the root's boundaries.
		tx0 = (min.x - origin.x) * invDirX;
		tx1 = (max.x - origin.x) * invDirX;
		ty0 = (min.y - origin.y) * invDirY;
		ty1 = (max.y - origin.y) * invDirY;
		tz0 = (min.z - origin.z) * invDirZ;
		tz1 = (max.z - origin.z) * invDirZ;

		// Check if the ray hits the octree.
		if(Math.max(Math.max(tx0, ty0), tz0) < Math.min(Math.min(tx1, ty1), tz1)) {

			// Find the intersecting octants.
			raycastOctant(octree.root, tx0, ty0, tz0, tx1, ty1, tz1, raycaster, intersects);

		}

	};

	/**
	 * A 3D box.
	 *
	 * @type {Box3}
	 * @private
	 */

	var b$3 = new Box3();

	/**
	 * Recursively calculates the depth of the given octree.
	 *
	 * @private
	 * @param {Octant} octant - An octant.
	 * @return {Number} The depth.
	 */

	function getDepth(octant) {

		var children = octant.children;

		var result = 0;
		var i, l, d;

		if(children !== null) {

			for(i = 0, l = children.length; i < l; ++i) {

				d = 1 + getDepth(children[i]);

				if(d > result) {

					result = d;

				}

			}

		}

		return result;

	}

	/**
	 * Recursively collects octants that lie inside the specified region.
	 *
	 * @private
	 * @param {Octant} octant - An octant.
	 * @param {Frustum|Box3} region - A region.
	 * @param {Octant[]} result - A list to be filled with octants that intersect with the region.
	 */

	function cull(octant, region, result) {

		var children = octant.children;

		var i, l;

		b$3.min = octant.min;
		b$3.max = octant.max;

		if(region.intersectsBox(b$3)) {

			if(children !== null) {

				for(i = 0, l = children.length; i < l; ++i) {

					cull(children[i], region, result);

				}

			} else {

				result.push(octant);

			}

		}

	}

	/**
	 * Recursively fetches all octants with the specified depth level.
	 *
	 * @private
	 * @param {Octant} octant - An octant.
	 * @param {Number} level - The target depth level.
	 * @param {Number} depth - The current depth level.
	 * @param {Octant[]} result - A list to be filled with the identified octants.
	 */

	function findOctantsByLevel(octant, level, depth, result) {

		var children = octant.children;

		var i, l;

		if(depth === level) {

			result.push(octant);

		} else if(children !== null) {

			++depth;

			for(i = 0, l = children.length; i < l; ++i) {

				findOctantsByLevel(children[i], level, depth, result);

			}

		}

	}

	/**
	 * An octree that subdivides space for fast spatial searches.
	 *
	 * @implements {Iterable}
	 */

	var Octree = function Octree(min, max) {

		/**
			 * The root octant.
			 *
			 * @type {Octant}
			 * @default null
			 */

		this.root = (min !== undefined && max !== undefined) ? new Octant(min, max) : null;

	};

	var prototypeAccessors$2 = { min: { configurable: true },max: { configurable: true },children: { configurable: true } };

	/**
		 * The lower bounds of the root octant.
		 *
		 * @type {Vector3}
		 */

	prototypeAccessors$2.min.get = function () { return this.root.min; };

	/**
		 * The upper bounds of the root octant.
		 *
		 * @type {Vector3}
		 */

	prototypeAccessors$2.max.get = function () { return this.root.max; };

	/**
		 * The children of the root octant.
		 *
		 * @type {Octant[]}
		 */

	prototypeAccessors$2.children.get = function () { return this.root.children; };

	/**
		 * Calculates the center of this octree.
		 *
		 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
		 * @return {Vector3} A vector that describes the center of this octree.
		 */

	Octree.prototype.getCenter = function getCenter (target) { return this.root.getCenter(target); };

	/**
		 * Calculates the size of this octree.
		 *
		 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
		 * @return {Vector3} A vector that describes the size of this octree.
		 */

	Octree.prototype.getDimensions = function getDimensions (target) { return this.root.getDimensions(target); };

	/**
		 * Calculates the current depth of this octree.
		 *
		 * @return {Number} The depth.
		 */

	Octree.prototype.getDepth = function getDepth$1 () { return getDepth(this.root); };

	/**
		 * Recursively collects octants that intersect with the specified region.
		 *
		 * @param {Frustum|Box3} region - A region.
		 * @return {Octant[]} The octants.
		 */

	Octree.prototype.cull = function cull$1 (region) {

		var result = [];

		cull(this.root, region, result);

		return result;

	};

	/**
		 * Fetches all octants with the specified depth level.
		 *
		 * @param {Number} level - The depth level.
		 * @return {Octant[]} The octants.
		 */

	Octree.prototype.findOctantsByLevel = function findOctantsByLevel$1 (level) {

		var result = [];

		findOctantsByLevel(this.root, level, 0, result);

		return result;

	};

	/**
		 * Finds the octants that intersect with the given ray. The intersecting
		 * octants are sorted by distance, closest first.
		 *
		 * @param {Raycaster} raycaster - A raycaster.
		 * @param {Octant[]} [intersects] - An optional target list to be filled with the intersecting octants.
		 * @return {Octant[]} The intersecting octants.
		 */

	Octree.prototype.raycast = function raycast (raycaster, intersects) {
			if ( intersects === void 0 ) intersects = [];


		OctreeRaycaster.intersectOctree(this, raycaster, intersects);

		return intersects;

	};

	/**
		 * Returns an iterator that traverses the octree and returns leaf nodes.
		 *
		 * When a cull region is provided, the iterator will only return leaves that
		 * intersect with that region.
		 *
		 * @param {Frustum|Box3} [region] - A cull region.
		 * @return {OctantIterator} An iterator.
		 */

	Octree.prototype.leaves = function leaves (region) {

		return new OctantIterator(this, region);

	};

	/**
		 * Returns an iterator that traverses the octree and returns all leaf nodes.
		 *
		 * @return {OctantIterator} An iterator.
		 */

	Octree.prototype[Symbol.iterator] = function () {

		return new OctantIterator(this);

	};

	Object.defineProperties( Octree.prototype, prototypeAccessors$2 );

	/**
	 * Core components.
	 *
	 * @module sparse-octree/core
	 */

	/**
	 * A point.
	 *
	 * @type {Vector3}
	 * @private
	 */

	var p = new Vector3();

	/**
	 * An octant that maintains points.
	 */

	var PointOctant = (function (Octant$$1) {
		function PointOctant(min, max) {

			Octant$$1.call(this, min, max);

			/**
			 * The points.
			 *
			 * @type {Vector3[]}
			 */

			this.points = null;

			/**
			 * Point data.
			 *
			 * @type {Array}
			 */

			this.data = null;

		}

		if ( Octant$$1 ) PointOctant.__proto__ = Octant$$1;
		PointOctant.prototype = Object.create( Octant$$1 && Octant$$1.prototype );
		PointOctant.prototype.constructor = PointOctant;

		/**
		 * Computes the distance squared from this octant to the given point.
		 *
		 * @param {Vector3} point - A point.
		 * @return {Number} The distance squared.
		 */

		PointOctant.prototype.distanceToSquared = function distanceToSquared (point) {

			var clampedPoint = p.copy(point).clamp(this.min, this.max);

			return clampedPoint.sub(point).lengthSquared();

		};

		/**
		 * Computes the distance squared from the center of this octant to the given
		 * point.
		 *
		 * @param {Vector3} point - A point.
		 * @return {Number} The distance squared.
		 */

		PointOctant.prototype.distanceToCenterSquared = function distanceToCenterSquared (point) {

			var center = this.getCenter(p);

			var dx = point.x - center.x;
			var dy = point.y - center.x;
			var dz = point.z - center.z;

			return dx * dx + dy * dy + dz * dz;

		};

		/**
		 * Checks if the given point lies inside this octant's boundaries.
		 *
		 * This method can also be used to check if this octant intersects a sphere by
		 * providing a radius as bias.
		 *
		 * @param {Vector3} point - A point.
		 * @param {Number} bias - A padding that extends the boundaries temporarily.
		 * @return {Boolean} Whether the given point lies inside this octant.
		 */

		PointOctant.prototype.contains = function contains (point, bias) {

			var min = this.min;
			var max = this.max;

			return (
				point.x >= min.x - bias &&
				point.y >= min.y - bias &&
				point.z >= min.z - bias &&
				point.x <= max.x + bias &&
				point.y <= max.y + bias &&
				point.z <= max.z + bias
			);

		};

		/**
		 * Redistributes existing points to child octants.
		 *
		 * @param {Number} bias - A proximity threshold.
		 */

		PointOctant.prototype.redistribute = function redistribute (bias) {

			var children = this.children;
			var points = this.points;
			var data = this.data;

			var i, j, il, jl;
			var child, point, entry;

			if(children !== null) {

				for(i = 0, il = points.length; i < il; ++i) {

					point = points[i];
					entry = data[i];

					for(j = 0, jl = children.length; j < jl; ++j) {

						child = children[j];

						if(child.contains(point, bias)) {

							if(child.points === null) {

								child.points = [];
								child.data = [];

							}

							child.points.push(point);
							child.data.push(entry);

							break;

						}

					}

				}

			}

			this.points = null;
			this.data = null;

		};

		/**
		 * Gathers all points from the children. The children are expected to be leaf
		 * octants and will be dropped afterwards.
		 */

		PointOctant.prototype.merge = function merge () {
			var this$1 = this;


			var children = this.children;

			var i, l;
			var child;

			if(children !== null) {

				this.points = [];
				this.data = [];

				for(i = 0, l = children.length; i < l; ++i) {

					child = children[i];

					if(child.points !== null) {

						(ref = this$1.points).push.apply(ref, child.points);
						(ref$1 = this$1.data).push.apply(ref$1, child.data);

					}

				}

				this.children = null;

			}
			var ref;
			var ref$1;

		};

		return PointOctant;
	}(Octant));

	/**
	 * A collection of ray-point intersection data.
	 */

	var RayPointIntersection = function RayPointIntersection(distance, distanceToRay, point, object) {
		if ( object === void 0 ) object = null;


		/**
			 * The distance from the origin of the ray to the point.
			 *
			 * @type {Number}
			 */

		this.distance = distance;

		/**
			 * The shortest distance from the point to the ray.
			 *
			 * @type {Number}
			 */

		this.distanceToRay = distanceToRay;

		/**
			 * The point.
			 *
			 * @type {Vector3}
			 */

		this.point = point;

		/**
			 * The point's data.
			 *
			 * @type {Object}
			 * @default null
			 */

		this.object = object;

	};

	/**
	 * A threshold for distance comparisons.
	 *
	 * @type {Number}
	 * @private
	 */

	var THRESHOLD = 1e-6;

	/**
	 * Recursively counts how many points are in the given octant.
	 *
	 * @private
	 * @param {Octant} octant - An octant.
	 * @return {Number} The amount of points.
	 */

	function countPoints(octant) {

		var children = octant.children;

		var result = 0;
		var i, l;

		if(children !== null) {

			for(i = 0, l = children.length; i < l; ++i) {

				result += countPoints(children[i]);

			}

		} else if(octant.points !== null) {

			result = octant.points.length;

		}

		return result;

	}

	/**
	 * Recursively places a point into the octree.
	 *
	 * @private
	 * @param {Vector3} point - A point.
	 * @param {Object} data - An object that the point represents.
	 * @param {Octree} octree - The octree.
	 * @param {Octant} octant - The current octant.
	 * @param {Number} depth - The current depth.
	 * @return {Boolean} Whether the operation was successful.
	 */

	function put(point, data, octree, octant, depth) {

		var children = octant.children;
		var exists = false;
		var done = false;
		var i, l;

		if(octant.contains(point, octree.bias)) {

			if(children === null) {

				if(octant.points === null) {

					octant.points = [];
					octant.data = [];

				} else {

					for(i = 0, l = octant.points.length; !exists && i < l; ++i) {

						exists = octant.points[i].equals(point);

					}

				}

				if(exists) {

					octant.data[i - 1] = data;
					done = true;

				} else if(octant.points.length < octree.maxPoints || depth === octree.maxDepth) {

					octant.points.push(point.clone());
					octant.data.push(data);
					++octree.pointCount;
					done = true;

				} else {

					octant.split();
					octant.redistribute(octree.bias);
					children = octant.children;

				}

			}

			if(children !== null) {

				++depth;

				for(i = 0, l = children.length; !done && i < l; ++i) {

					done = put(point, data, octree, children[i], depth);

				}

			}

		}

		return done;

	}

	/**
	 * Recursively finds a point in the octree and removes it.
	 *
	 * @private
	 * @param {Vector3} point - A point.
	 * @param {Octree} octree - The octree.
	 * @param {Octant} octant - The current octant.
	 * @param {Octant} parent - The parent of the current octant.
	 * @return {Object} The data entry of the removed point or null if it didn't exist.
	 */

	function remove(point, octree, octant, parent) {

		var children = octant.children;

		var result = null;

		var i, l;
		var points, data, last;

		if(octant.contains(point, octree.bias)) {

			if(children !== null) {

				for(i = 0, l = children.length; result === null && i < l; ++i) {

					result = remove(point, octree, children[i], octant);

				}

			} else if(octant.points !== null) {

				points = octant.points;
				data = octant.data;

				for(i = 0, l = points.length; i < l; ++i) {

					if(points[i].equals(point)) {

						last = l - 1;
						result = data[i];

						// If the point is NOT the last one in the array:
						if(i < last) {

							// Overwrite with the last point and data entry.
							points[i] = points[last];
							data[i] = data[last];

						}

						// Drop the last entry.
						points.pop();
						data.pop();

						--octree.pointCount;

						if(parent !== null && countPoints(parent) <= octree.maxPoints) {

							parent.merge();

						}

						break;

					}

				}

			}

		}

		return result;

	}

	/**
	 * Recursively finds a point in the octree and fetches the associated data.
	 *
	 * @private
	 * @param {Vector3} point - A point.
	 * @param {Octree} octree - The octree.
	 * @param {Octant} octant - The current octant octant.
	 * @return {Object} The data entry that is associated with the given point or null if it doesn't exist.
	 */

	function fetch(point, octree, octant) {

		var children = octant.children;

		var result = null;

		var i, l;
		var points;

		if(octant.contains(point, octree.bias)) {

			if(children !== null) {

				for(i = 0, l = children.length; result === null && i < l; ++i) {

					result = fetch(point, octree, children[i]);

				}

			} else {

				points = octant.points;

				for(i = 0, l = points.length; result === null && i < l; ++i) {

					if(point.distanceToSquared(points[i]) <= THRESHOLD) {

						result = octant.data[i];

					}

				}

			}

		}

		return result;

	}

	/**
	 * Recursively moves an existing point to a new position.
	 *
	 * @private
	 * @param {Vector3} point - The point.
	 * @param {Vector3} position - The new position.
	 * @param {Octree} octree - The octree.
	 * @param {Octant} octant - The current octant.
	 * @param {Octant} octant - The parent of the current octant.
	 * @return {Object} The data entry of the updated point or null if it didn't exist.
	 */

	function move(point, position, octree, octant, parent, depth) {

		var children = octant.children;

		var result = null;

		var i, l;
		var points;

		if(octant.contains(point, octree.bias)) {

			if(octant.contains(position, octree.bias)) {

				// The point and the new position both fall into the current octant.
				if(children !== null) {

					++depth;

					for(i = 0, l = children.length; result === null && i < l; ++i) {

						result = move(point, position, octree, children[i], octant, depth);

					}

				} else {

					// No divergence - the point can be updated in place.
					points = octant.points;

					for(i = 0, l = points.length; i < l; ++i) {

						if(point.distanceToSquared(points[i]) <= THRESHOLD) {

							// The point exists! Update its position.
							points[i].copy(position);
							result = octant.data[i];

							break;

						}

					}

				}

			} else {

				// Retrieve the point and remove it.
				result = remove(point, octree, octant, parent);

				// Go back to the parent octant and add the updated point.
				put(position, result, octree, parent, depth - 1);

			}

		}

		return result;

	}

	/**
	 * Recursively finds the closest point to the given one.
	 *
	 * @private
	 * @param {Vector3} point - The point.
	 * @param {Number} maxDistance - The maximum distance.
	 * @param {Boolean} skipSelf - Whether a point that is exactly at the given position should be skipped.
	 * @param {Octant} octant - The current octant.
	 * @return {Object} An object representing the nearest point or null if there is none. The object has a point and a data property.
	 */

	function findNearestPoint(point, maxDistance, skipSelf, octant) {

		var points = octant.points;
		var children = octant.children;

		var result = null;
		var bestDist = maxDistance;

		var i, l;
		var p, distSq;

		var sortedChildren;
		var child, childResult;

		if(children !== null) {

			// Sort the children.
			sortedChildren = children.map(function(child) {

				// Precompute distances.
				return {
					octant: child,
					distance: child.distanceToCenterSquared(point)
				};

			}).sort(function(a, b) {

				// Smallest distance to the point first, ASC.
				return a.distance - b.distance;

			});

			// Traverse from closest to furthest.
			for(i = 0, l = sortedChildren.length; i < l; ++i) {

				// Unpack octant.
				child = sortedChildren[i].octant;

				if(child.contains(point, bestDist)) {

					childResult = findNearestPoint(point, bestDist, skipSelf, child);

					if(childResult !== null) {

						distSq = childResult.point.distanceToSquared(point);

						if((!skipSelf || distSq > 0.0) && distSq < bestDist) {

							bestDist = distSq;
							result = childResult;

						}

					}

				}

			}

		} else if(points !== null) {

			for(i = 0, l = points.length; i < l; ++i) {

				p = points[i];
				distSq = point.distanceToSquared(p);

				if((!skipSelf || distSq > 0.0) && distSq < bestDist) {

					bestDist = distSq;

					result = {
						point: p.clone(),
						data: octant.data[i]
					};

				}

			}

		}

		return result;

	}

	/**
	 * Recursively finds points that are inside the specified radius around a given
	 * position.
	 *
	 * @private
	 * @param {Vector3} point - A position.
	 * @param {Number} radius - A radius.
	 * @param {Boolean} skipSelf - Whether a point that is exactly at the given position should be skipped.
	 * @param {Octant} octant - The current octant.
	 * @param {Array} result - An array to be filled with objects, each containing a point and a data property.
	 */

	function findPoints(point, radius, skipSelf, octant, result) {

		var points = octant.points;
		var children = octant.children;
		var rSq = radius * radius;

		var i, l;

		var p, distSq;
		var child;

		if(children !== null) {

			for(i = 0, l = children.length; i < l; ++i) {

				child = children[i];

				if(child.contains(point, radius)) {

					findPoints(point, radius, skipSelf, child, result);

				}

			}

		} else if(points !== null) {

			for(i = 0, l = points.length; i < l; ++i) {

				p = points[i];
				distSq = point.distanceToSquared(p);

				if((!skipSelf || distSq > 0.0) && distSq <= rSq) {

					result.push({
						point: p.clone(),
						data: octant.data[i]
					});

				}

			}

		}

	}

	/**
	 * An octree that manages points.
	 */

	var PointOctree = (function (Octree$$1) {
		function PointOctree(min, max, bias, maxPoints, maxDepth) {
			if ( bias === void 0 ) bias = 0.0;
			if ( maxPoints === void 0 ) maxPoints = 8;
			if ( maxDepth === void 0 ) maxDepth = 8;


			Octree$$1.call(this);

			/**
			 * The root octant.
			 *
			 * @type {PointOctant}
			 */

			this.root = new PointOctant(min, max);

			/**
			 * An octant boundary bias.
			 *
			 * The octree is considered "loose" with a bias greater than 0.
			 *
			 * @type {Number}
			 * @private
			 * @default 0.0
			 */

			this.bias = Math.max(0.0, bias);

			/**
			 * Number of points per octant before a split occurs.
			 *
			 * This value works together with the maximum depth as a secondary limiting
			 * factor. Smaller values cause splits to occur earlier which results in a
			 * faster and deeper tree growth.
			 *
			 * @type {Number}
			 * @private
			 * @default 8
			 */

			this.maxPoints = Math.max(1, Math.round(maxPoints));

			/**
			 * The maximum tree depth level.
			 *
			 * It's possible to use Infinity, but keep in mind that allowing infinitely
			 * small octants can have a severely negative impact on performance.
			 * Finding a value that works best for a specific scene is advisable.
			 *
			 * @type {Number}
			 * @private
			 * @default 8
			 */

			this.maxDepth = Math.max(0, Math.round(maxDepth));

			/**
			 * The amount of points that are currently in this octree.
			 *
			 * @type {Number}
			 */

			this.pointCount = 0;

		}

		if ( Octree$$1 ) PointOctree.__proto__ = Octree$$1;
		PointOctree.prototype = Object.create( Octree$$1 && Octree$$1.prototype );
		PointOctree.prototype.constructor = PointOctree;

		/**
		 * Counts how many points are in the given octant.
		 *
		 * @param {Octant} octant - An octant.
		 * @return {Number} The amount of points.
		 */

		PointOctree.prototype.countPoints = function countPoints$1 (octant) {

			return countPoints(octant);

		};

		/**
		 * Puts a point into the octree.
		 *
		 * @param {Vector3} point - A point. If it's already in the octree, the data entry will be updated.
		 * @param {Object} data - A data object that belongs to the point.
		 * @return {Boolean} Whether the operation was successful.
		 */

		PointOctree.prototype.put = function put$1 (point, data) {

			return put(point, data, this, this.root, 0);

		};

		/**
		 * Removes a point from the tree.
		 *
		 * @param {Vector3} point - A point.
		 * @return {Object} The data entry of the removed point or null if it didn't exist.
		 */

		PointOctree.prototype.remove = function remove$1 (point) {

			return remove(point, this, this.root, null);

		};

		/**
		 * Retrieves the data of the specified point.
		 *
		 * @param {Vector3} point - A position.
		 * @return {Object} The data entry that is associated with the given point or null if it doesn't exist.
		 */

		PointOctree.prototype.fetch = function fetch$1 (point) {

			return fetch(point, this, this.root);

		};

		/**
		 * Moves an existing point to a new position. Has no effect if the point
		 * doesn't exist.
		 *
		 * @param {Vector3} point - The point.
		 * @param {Vector3} position - The new position.
		 * @return {Object} The data entry of the updated point or null if it didn't exist.
		 */

		PointOctree.prototype.move = function move$1 (point, position) {

			return move(point, position, this, this.root, null, 0);

		};

		/**
		 * Finds the closest point to the given one.
		 *
		 * @param {Vector3} point - A point.
		 * @param {Number} [maxDistance=Infinity] - An upper limit for the distance between the points.
		 * @param {Boolean} [skipSelf=false] - Whether a point that is exactly at the given position should be skipped.
		 * @return {Object} An object representing the nearest point or null if there is none. The object has a point and a data property.
		 */

		PointOctree.prototype.findNearestPoint = function findNearestPoint$1 (point, maxDistance, skipSelf) {
			if ( maxDistance === void 0 ) maxDistance = Infinity;
			if ( skipSelf === void 0 ) skipSelf = false;


			return findNearestPoint(point, maxDistance, skipSelf, this.root);

		};

		/**
		 * Finds points that are in the specified radius around the given position.
		 *
		 * @param {Vector3} point - A position.
		 * @param {Number} radius - A radius.
		 * @param {Boolean} [skipSelf=false] - Whether a point that is exactly at the given position should be skipped.
		 * @return {Array} An array of objects, each containing a point and a data property.
		 */

		PointOctree.prototype.findPoints = function findPoints$1 (point, radius, skipSelf) {
			if ( skipSelf === void 0 ) skipSelf = false;


			var result = [];

			findPoints(point, radius, skipSelf, this.root, result);

			return result;

		};

		/**
		 * Finds the points that intersect with the given ray.
		 *
		 * @param {Raycaster} raycaster - The raycaster.
		 * @param {Array} [intersects] - An array to be filled with the intersecting points.
		 * @return {RayPointIntersection[]} The intersecting points.
		 */

		PointOctree.prototype.raycast = function raycast (raycaster, intersects) {
			if ( intersects === void 0 ) intersects = [];


			var octants = Octree$$1.prototype.raycast.call(this, raycaster);

			if(octants.length > 0) {

				// Collect intersecting points.
				this.testPoints(octants, raycaster, intersects);

			}

			return intersects;

		};

		/**
		 * Collects points that intersect with the given ray.
		 *
		 * @param {Octant[]} octants - An array containing octants that intersect with the ray.
		 * @param {Raycaster} raycaster - The raycaster.
		 * @param {Array} intersects - An array to be filled with intersecting points.
		 */

		PointOctree.prototype.testPoints = function testPoints (octants, raycaster, intersects) {

			var threshold = raycaster.params.Points.threshold;
			var thresholdSq = threshold * threshold;

			var intersectPoint;
			var distance, distanceToRay;
			var rayPointDistanceSq;

			var i, j, il, jl;
			var octant, points, point;

			for(i = 0, il = octants.length; i < il; ++i) {

				octant = octants[i];
				points = octant.points;

				if(points !== null) {

					for(j = 0, jl = points.length; j < jl; ++j) {

						point = points[j];
						rayPointDistanceSq = raycaster.ray.distanceSqToPoint(point);

						if(rayPointDistanceSq < thresholdSq) {

							intersectPoint = raycaster.ray.closestPointToPoint(point);
							distance = raycaster.ray.origin.distanceTo(intersectPoint);

							if(distance >= raycaster.near && distance <= raycaster.far) {

								distanceToRay = Math.sqrt(rayPointDistanceSq);

								intersects.push(new RayPointIntersection(
									distance,
									distanceToRay,
									intersectPoint,
									octant.data[j]
								));

							}

						}

					}

				}

			}

		};

		return PointOctree;
	}(Octree));

	/**
	 * Point-oriented octree components.
	 *
	 * @module sparse-octree/points
	 */

	/**
	 * A collection of octree utility functions.
	 */

	var OctreeUtils = function OctreeUtils () {};

	OctreeUtils.recycleOctants = function recycleOctants (octant, octants) {

		var a = new Vector3();
		var b = new Vector3();
		var c = new Vector3();

		var min = octant.min;
		var mid = octant.getCenter();
		var halfDimensions = octant.getDimensions().multiplyScalar(0.5);

		var children = octant.children;
		var l = octants.length;

		var i, j;
		var combination, candidate;

		for(i = 0; i < 8; ++i) {

			combination = pattern[i];

			b.addVectors(min, a.fromArray(combination).multiply(halfDimensions));
			c.addVectors(mid, a.fromArray(combination).multiply(halfDimensions));

			// Find an octant that matches the current combination.
			for(j = 0; j < l; ++j) {

				candidate = octants[j];

				if(candidate !== null && b.equals(candidate.min) && c.equals(candidate.max)) {

					children[i] = candidate;
					octants[j] = null;

					break;

				}

			}

		}

	};

	/**
	 * Octree utilities.
	 *
	 * @module sparse-octree/utils
	 */

	/**
	 * Exposure of the library components.
	 *
	 * @module sparse-octree
	 */

	exports.CubicOctant = CubicOctant;
	exports.edges = edges;
	exports.Octant = Octant;
	exports.Octree = Octree;
	exports.OctantIterator = OctantIterator;
	exports.OctreeRaycaster = OctreeRaycaster;
	exports.pattern = pattern;
	exports.PointOctant = PointOctant;
	exports.PointOctree = PointOctree;
	exports.RayPointIntersection = RayPointIntersection;
	exports.OctreeUtils = OctreeUtils;

	Object.defineProperty(exports, '__esModule', { value: true });

})));
