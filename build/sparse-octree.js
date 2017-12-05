/**
 * sparse-octree v4.1.1 build Dec 05 2017
 * https://github.com/vanruesc/sparse-octree
 * Copyright 2017 Raoul van Rüschen, Zlib
 */

(function (global, factory) {
	typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports, require('three')) :
	typeof define === 'function' && define.amd ? define(['exports', 'three'], factory) :
	(factory((global.SPARSEOCTREE = {}),global.three));
}(this, (function (exports,three) { 'use strict';

	/**
	 * A vector.
	 *
	 * @type {Vector3}
	 * @private
	 */

	var c$1 = new three.Vector3();

	/**
	 * An octant.
	 */

	var Octant = function Octant(min, max) {
		if ( min === void 0 ) min = new three.Vector3();
		if ( max === void 0 ) max = new three.Vector3();


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
			if ( target === void 0 ) target = new three.Vector3();


		return target.addVectors(this.min, this.max).multiplyScalar(0.5);

	};

	/**
		 * Computes the size of this octant.
		 *
		 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
		 * @return {Vector3} A vector that describes the size of this octant.
		 */

	Octant.prototype.getDimensions = function getDimensions (target) {
			if ( target === void 0 ) target = new three.Vector3();


		return target.subVectors(this.max, this.min);

	};

	/**
		 * Splits this octant into eight smaller ones.
		 */

	Octant.prototype.split = function split () {
			var this$1 = this;


		var min = this.min;
		var max = this.max;
		var mid = this.getCenter(c$1);

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

				new three.Vector3(
					(combination[0] === 0) ? min.x : mid.x,
					(combination[1] === 0) ? min.y : mid.y,
					(combination[2] === 0) ? min.z : mid.z
				),

				new three.Vector3(
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

	var c = new three.Vector3();

	/**
	 * A cubic octant.
	 */

	var CubicOctant = function CubicOctant(min, size) {
		if ( min === void 0 ) min = new three.Vector3();
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

	prototypeAccessors.max.get = function () {

		return this.min.clone().addScalar(this.size);

	};

	/**
		 * Computes the center of this octant.
		 *
		 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
		 * @return {Vector3} A vector that describes the center of this octant.
		 */

	CubicOctant.prototype.getCenter = function getCenter (target) {
			if ( target === void 0 ) target = new three.Vector3();


		return target.copy(this.min).addScalar(this.size * 0.5);

	};

	/**
		 * Returns the size of this octant as a vector.
		 *
		 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
		 * @return {Vector3} A vector that describes the size of this octant.
		 */

	CubicOctant.prototype.getDimensions = function getDimensions (target) {
			if ( target === void 0 ) target = new three.Vector3();


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

				new three.Vector3(
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

	var b$1 = new three.Box3();

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

			b$1.min = root.min;
			b$1.max = root.max;

			if(!this.cull || this.region.intersectsBox(b$1)) {

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

						b$1.min = child.min;
						b$1.max = child.max;

						if(!region.intersectsBox(b$1)) {

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
	 * A list of vectors.
	 *
	 * @type {Vector3[]}
	 * @private
	 * @final
	 */

	var v = [
		new three.Vector3(),
		new three.Vector3(),
		new three.Vector3()
	];

	/**
	 * A box.
	 *
	 * @type {Box3}
	 * @private
	 * @final
	 */

	var b$2 = new three.Box3();

	/**
	 * A ray.
	 *
	 * @type {Ray}
	 * @private
	 * @final
	 */

	var r = new three.Ray();

	/**
	 * A lookup-table containing octant ids. Used to determine the exit plane from
	 * an octant.
	 *
	 * @type {Uint8Array[]}
	 * @private
	 * @final
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
	 * A byte that stores raycasting flags.
	 *
	 * @type {Number}
	 * @private
	 */

	var flags = 0;

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
			if(tym < tx0) {

				entry |= 2;

			}

			if(tzm < tx0) {

				entry |= 1;

			}

		} else if(ty0 > tz0) {

			// XZ-plane.
			if(txm < ty0) {

				entry |= 4;

			}

			if(tzm < ty0) {

				entry |= 1;

			}

		} else {

			// XY-plane.
			if(txm < tz0) {

				entry |= 4;

			}

			if(tym < tz0) {

				entry |= 2;

			}

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
							raycastOctant(children[flags], tx0, ty0, tz0, txm, tym, tzm, raycaster, intersects);
							currentOctant = findNextOctant(currentOctant, txm, tym, tzm);
							break;

						case 1:
							raycastOctant(children[flags ^ 1], tx0, ty0, tzm, txm, tym, tz1, raycaster, intersects);
							currentOctant = findNextOctant(currentOctant, txm, tym, tz1);
							break;

						case 2:
							raycastOctant(children[flags ^ 2], tx0, tym, tz0, txm, ty1, tzm, raycaster, intersects);
							currentOctant = findNextOctant(currentOctant, txm, ty1, tzm);
							break;

						case 3:
							raycastOctant(children[flags ^ 3], tx0, tym, tzm, txm, ty1, tz1, raycaster, intersects);
							currentOctant = findNextOctant(currentOctant, txm, ty1, tz1);
							break;

						case 4:
							raycastOctant(children[flags ^ 4], txm, ty0, tz0, tx1, tym, tzm, raycaster, intersects);
							currentOctant = findNextOctant(currentOctant, tx1, tym, tzm);
							break;

						case 5:
							raycastOctant(children[flags ^ 5], txm, ty0, tzm, tx1, tym, tz1, raycaster, intersects);
							currentOctant = findNextOctant(currentOctant, tx1, tym, tz1);
							break;

						case 6:
							raycastOctant(children[flags ^ 6], txm, tym, tz0, tx1, ty1, tzm, raycaster, intersects);
							currentOctant = findNextOctant(currentOctant, tx1, ty1, tzm);
							break;

						case 7:
							raycastOctant(children[flags ^ 7], txm, tym, tzm, tx1, ty1, tz1, raycaster, intersects);
							// Far top right octant. No other octants can be reached from here.
							currentOctant = 8;
							break;

					}

				} while(currentOctant < 8);

			}

		}

	}

	/**
	 * An octree raycaster.
	 *
	 * Based on:
	 *  "An Efficient Parametric Algorithm for Octree Traversal"
	 *  by J. Revelles et al. (2000).
	 */

	var OctreeRaycaster = function OctreeRaycaster () {};

	OctreeRaycaster.intersectOctree = function intersectOctree (octree, raycaster, intersects) {

		// Translate the octree extents to the scene origin.
		var min = b$2.min.set(0, 0, 0);
		var max = b$2.max.subVectors(octree.max, octree.min);

		var dimensions = octree.getDimensions(v[0]);
		var halfDimensions = v[1].copy(dimensions).multiplyScalar(0.5);

		var origin = r.origin.copy(raycaster.ray.origin);
		var direction = r.direction.copy(raycaster.ray.direction);

		var invDirX, invDirY, invDirZ;
		var tx0, tx1, ty0, ty1, tz0, tz1;

		// Translate the ray to the center of the octree.
		origin.sub(octree.getCenter(v[2])).add(halfDimensions);

		// Reset all flags.
		flags = 0;

		// Handle rays with negative directions.
		if(direction.x < 0.0) {

			origin.x = dimensions.x - origin.x;
			direction.x = -direction.x;
			flags |= 4;

		}

		if(direction.y < 0.0) {

			origin.y = dimensions.y - origin.y;
			direction.y = -direction.y;
			flags |= 2;

		}

		if(direction.z < 0.0) {

			origin.z = dimensions.z - origin.z;
			direction.z = -direction.z;
			flags |= 1;

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

	var b = new three.Box3();

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

		b.min = octant.min;
		b.max = octant.max;

		if(region.intersectsBox(b)) {

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

	var prototypeAccessors$1 = { min: { configurable: true },max: { configurable: true },children: { configurable: true } };

	/**
		 * The lower bounds of the root octant.
		 *
		 * @type {Vector3}
		 */

	prototypeAccessors$1.min.get = function () {

		return this.root.min;

	};

	/**
		 * The upper bounds of the root octant.
		 *
		 * @type {Vector3}
		 */

	prototypeAccessors$1.max.get = function () {

		return this.root.max;

	};

	/**
		 * The children of the root octant.
		 *
		 * @type {Octant[]}
		 */

	prototypeAccessors$1.children.get = function () {

		return this.root.children;

	};

	/**
		 * Calculates the center of this octree.
		 *
		 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
		 * @return {Vector3} A vector that describes the center of this octree.
		 */

	Octree.prototype.getCenter = function getCenter (target) {

		return this.root.getCenter(target);

	};

	/**
		 * Calculates the size of this octree.
		 *
		 * @param {Vector3} [target] - A target vector. If none is provided, a new one will be created.
		 * @return {Vector3} A vector that describes the size of this octree.
		 */

	Octree.prototype.getDimensions = function getDimensions (target) {

		return this.root.getDimensions(target);

	};

	/**
		 * Calculates the current depth of this octree.
		 *
		 * @return {Number} The depth.
		 */

	Octree.prototype.getDepth = function getDepth$1 () {

		return getDepth(this.root);

	};

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

	Object.defineProperties( Octree.prototype, prototypeAccessors$1 );

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

	var p = new three.Vector3();

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

			return clampedPoint.sub(point).lengthSq();

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
	 * A box.
	 *
	 * @type {Box3}
	 * @private
	 * @final
	 */

	var b$3 = new three.Box3();

	/**
	 * A vector.
	 *
	 * @type {Vector3}
	 * @private
	 * @final
	 */

	var c$2 = new three.Vector3();

	/**
	 * A vector.
	 *
	 * @type {Vector3}
	 * @private
	 * @final
	 */

	var u = new three.Vector3();

	/**
	 * A vector.
	 *
	 * @type {Vector3}
	 * @private
	 * @final
	 */

	var v$1 = new three.Vector3();

	/**
	 * A collection of octree utility functions.
	 */

	var OctreeUtils = function OctreeUtils () {};

	OctreeUtils.recycleOctants = function recycleOctants (octant, octants) {

		var min = octant.min;
		var mid = octant.getCenter(u);
		var halfDimensions = octant.getDimensions(v$1).multiplyScalar(0.5);

		var children = octant.children;
		var l = octants.length;

		var i, j;
		var combination, candidate;

		for(i = 0; i < 8; ++i) {

			combination = pattern[i];

			b$3.min.addVectors(min, c$2.fromArray(combination).multiply(halfDimensions));
			b$3.max.addVectors(mid, c$2.fromArray(combination).multiply(halfDimensions));

			// Find an octant that matches the current combination.
			for(j = 0; j < l; ++j) {

				candidate = octants[j];

				if(candidate !== null && b$3.min.equals(candidate.min) && b$3.max.equals(candidate.max)) {

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
