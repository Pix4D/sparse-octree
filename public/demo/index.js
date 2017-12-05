(function (three,dat,Stats) {
	'use strict';

	dat = dat && dat.hasOwnProperty('default') ? dat['default'] : dat;
	Stats = Stats && Stats.hasOwnProperty('default') ? Stats['default'] : Stats;

	/**
	 * An octree helper.
	 */

	var OctreeHelper = (function (Group$$1) {
		function OctreeHelper(octree) {
			if ( octree === void 0 ) octree = null;


			Group$$1.call(this);

			/**
			 * The name of this object.
			 */

			this.name = "OctreeHelper";

			/**
			 * The octree.
			 *
			 * @type {Octree}
			 * @default null
			 */

			this.octree = octree;

			this.update();

		}

		if ( Group$$1 ) OctreeHelper.__proto__ = Group$$1;
		OctreeHelper.prototype = Object.create( Group$$1 && Group$$1.prototype );
		OctreeHelper.prototype.constructor = OctreeHelper;

		/**
		 * Creates octant geometry.
		 *
		 * @private
		 * @param {Iterator} octants - An octant iterator.
		 * @param {Number} octantCount - The size of the given sequence.
		 */

		OctreeHelper.prototype.createLineSegments = function createLineSegments (octants, octantCount) {

			var maxOctants = (Math.pow(2, 16) / 8) - 1;
			var group = new Group$$1();

			var material = new three.LineBasicMaterial({
				color: 0xffffff * Math.random()
			});

			var result;
			var vertexCount;
			var length;

			var indices, positions;
			var octant, min, max;
			var geometry;

			var i, j, c, d, n;
			var corner, edge;

			// Create geometry in multiple runs to limit the amount of vertices.
			for(i = 0, length = 0, n = Math.ceil(octantCount / maxOctants); n > 0; --n) {

				length += (octantCount < maxOctants) ? octantCount : maxOctants;
				octantCount -= maxOctants;

				vertexCount = length * 8;
				indices = new Uint16Array(vertexCount * 3);
				positions = new Float32Array(vertexCount * 3);

				// Continue where the previous run left off.
				for(c = 0, d = 0, result = octants.next(); !result.done && i < length;) {

					octant = result.value;
					min = octant.min;
					max = octant.max;

					// Create line connections based on the current vertex count.
					for(j = 0; j < 12; ++j) {

						edge = edges[j];

						indices[d++] = c + edge[0];
						indices[d++] = c + edge[1];

					}

					// Create the vertices.
					for(j = 0; j < 8; ++j, ++c) {

						corner = corners[j];

						positions[c * 3] = (corner[0] === 0) ? min.x : max.x;
						positions[c * 3 + 1] = (corner[1] === 0) ? min.y : max.y;
						positions[c * 3 + 2] = (corner[2] === 0) ? min.z : max.z;

					}

					if(++i < length) {

						result = octants.next();

					}

				}

				geometry = new three.BufferGeometry();
				geometry.setIndex(new three.BufferAttribute(indices, 1));
				geometry.addAttribute("position", new three.BufferAttribute(positions, 3));

				group.add(new three.LineSegments(geometry, material));

			}

			this.add(group);

		};

		/**
		 * Updates the helper geometry.
		 */

		OctreeHelper.prototype.update = function update () {
			var this$1 = this;


			var depth = (this.octree !== null) ? this.octree.getDepth() : -1;

			var level = 0;
			var result;

			// Remove existing geometry.
			this.dispose();

			while(level <= depth) {

				result = this$1.octree.findOctantsByLevel(level);

				this$1.createLineSegments(
					result[Symbol.iterator](),
					(typeof result.size === "number") ? result.size : result.length
				);

				++level;

			}

		};

		/**
		 * Destroys this helper.
		 */

		OctreeHelper.prototype.dispose = function dispose () {
			var this$1 = this;


			var groups = this.children;

			var group, children;
			var i, j, il, jl;

			for(i = 0, il = groups.length; i < il; ++i) {

				group = groups[i];
				children = group.children;

				for(j = 0, jl = children.length; j < jl; ++j) {

					children[j].geometry.dispose();
					children[j].material.dispose();

				}

				while(children.length > 0) {

					group.remove(children[0]);

				}

			}

			while(groups.length > 0) {

				this$1.remove(groups[0]);

			}

		};

		return OctreeHelper;
	}(three.Group));

	/**
	 * A binary pattern that describes the corners of an octant:
	 *
	 * ```text
	 *    3____7
	 *  2/___6/|
	 *  | 1__|_5
	 *  0/___4/
	 * ```
	 *
	 * @type {Uint8Array[]}
	 */

	var corners = [

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
	 * Exposure of the library components.
	 *
	 * @module octree-helper
	 */

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

	/**
	 * A demo base class.
	 */

	var Demo = function Demo(renderer) {

		/**
			 * A renderer.
			 *
			 * @type {WebGLRenderer}
			 */

		this.renderer = renderer;

		/**
			 * A loading manager.
			 *
			 * @type {LoadingManager}
			 */

		this.loadingManager = new three.LoadingManager();

		/**
			 * An asset map.
			 *
			 * @type {Map}
			 */

		this.assets = null;

		/**
			 * A scene.
			 *
			 * @type {Scene}
			 */

		this.scene = new three.Scene();
		this.scene.fog = new three.FogExp2(0x0d0d0d, 0.0025);

		/**
			 * A camera.
			 *
			 * @type {PerspectiveCamera}
			 */

		this.camera = new three.PerspectiveCamera(50, window.innerWidth / window.innerHeight, 0.1, 2000);

		/**
			 * Camera controls.
			 *
			 * @type {OrbitControls}
			 */

		this.controls = null;

	};

	/**
		 * Loads the demo. Override this method to load assets.
		 *
		 * @param {Function} callback - Call this function when all assets have been loaded.
		 */

	Demo.prototype.load = function load (callback) {

		callback();

	};

	/**
		 * Creates the scene.
		 */

	Demo.prototype.initialise = function initialise () {};

	/**
		 * Renders this demo.
		 *
		 * @param {Number} delta - The time since the last frame in seconds.
		 */

	Demo.prototype.render = function render (delta) {};

	/**
		 * Registers configuration options.
		 *
		 * @param {GUI} gui - A GUI.
		 */

	Demo.prototype.configure = function configure (gui) {};

	/**
		 * Resets this demo.
		 *
		 * @return {Demo} This demo.
		 */

	Demo.prototype.reset = function reset () {

		var fog = this.scene.fog;

		this.scene = new three.Scene();
		this.scene.fog = fog;

		if(this.controls !== null) {

			this.controls.dispose();
			this.controls = null;

		}

		return this;

	};

	/**
	 * A mouse position.
	 *
	 * @type {Vector2}
	 * @private
	 */

	var mouse = new three.Vector2();

	/**
	 * An octree raycaster.
	 *
	 * @implements {EventListener}
	 */

	var OctreeRaycaster$1 = (function (Raycaster$$1) {
		function OctreeRaycaster(octree, camera, object) {

			Raycaster$$1.call(this);

			/**
			 * A picking accuracy threshold for points.
			 */

			this.params.Points.threshold = 1e-1;

			/**
			 * An octree.
			 *
			 * @type {Octree}
			 * @private
			 */

			this.octree = octree;

			/**
			 * A camera.
			 *
			 * @type {PerspectiveCamera}
			 * @private
			 */

			this.camera = camera;

			/**
			 * An object to raycast with a brute force approach.
			 *
			 * @type {Object3D}
			 */

			this.object = object;

			/**
			 * Indicates whether the frustum culling is active.
			 *
			 * @type {Boolean}
			 * @default false
			 */

			this.enabled = true;

			/**
			 * A delta time.
			 *
			 * @type {String}
			 */

			this.delta = "";

			/**
			 * A selected object.
			 *
			 * @type {Object3D}
			 * @private
			 */

			this.selectedObject = null;

			/**
			 * The currently selected point.
			 *
			 * @type {Mesh}
			 */

			this.selectedPoint = new three.Mesh(
				new three.SphereBufferGeometry(0.2, 16, 16),
				new three.MeshBasicMaterial({
					transparent: true,
					color: 0x00ccff,
					opacity: 0.75
				})
			);

			this.selectedPoint.visible = false;

		}

		if ( Raycaster$$1 ) OctreeRaycaster.__proto__ = Raycaster$$1;
		OctreeRaycaster.prototype = Object.create( Raycaster$$1 && Raycaster$$1.prototype );
		OctreeRaycaster.prototype.constructor = OctreeRaycaster;

		/**
		 * Raycasts on mouse move events.
		 *
		 * @param {Event} event - A worker message event.
		 */

		OctreeRaycaster.prototype.handleEvent = function handleEvent (event) {

			switch(event.type) {

				case "mousemove":
					this.raycast(event);
					break;

			}

		};

		/**
		 * Raycasts the octree.
		 *
		 * @param {Event} event - An event.
		 */

		OctreeRaycaster.prototype.raycast = function raycast (event) {

			var intersects;
			var t0, t, x;

			mouse.x = (event.clientX / window.innerWidth) * 2 - 1;
			mouse.y = -(event.clientY / window.innerHeight) * 2 + 1;

			this.setFromCamera(mouse, this.camera);

			if(this.enabled) {

				// Use the octree raycasting capabilities.
				t0 = performance.now();
				intersects = this.octree.raycast(this);
				t = performance.now();

			} else {

				// Brute force alternative.
				t0 = performance.now();
				intersects = this.intersectObjects(this.object.children);
				t = performance.now();

			}

			this.delta = (t - t0).toFixed(2) + " ms";

			if(this.selectedObject !== null) {

				this.selectedObject.material.color.setHex(0xc00000);
				this.selectedObject = null;
				this.selectedPoint.visible = false;

			}

			if(intersects.length > 0) {

				x = intersects[0];

				if(x.object !== undefined) {

					this.selectedObject = x.object;
					this.selectedObject.material.color.setHex(0xccff00);

					this.selectedPoint.visible = x.object.parent.visible;
					this.selectedPoint.position.copy(x.point);

				} else {

					console.warn(intersects);

				}

			}

		};

		/**
		 * Registers configuration options.
		 *
		 * @param {GUI} gui - A GUI.
		 */

		OctreeRaycaster.prototype.configure = function configure (gui) {

			var folder = gui.addFolder("Raycasting");

			folder.add(this, "enabled");
			folder.add(this, "delta").listen();

			folder.open();

		};

		return OctreeRaycaster;
	}(three.Raycaster));

	/**
	 * A matrix.
	 *
	 * @type {Matrix4}
	 * @private
	 */

	var matrix4 = new three.Matrix4();

	/**
	 * A frustum.
	 *
	 * @type {Frustum}
	 * @private
	 */

	var frustum = new three.Frustum();

	/**
	 * A frustum-based octree culler.
	 */

	var FrustumCuller = function FrustumCuller(octree, scene) {

		/**
			 * An octree.
			 *
			 * @type {Octree}
			 * @private
			 */

		this.octree = octree;

		/**
			 * A scene.
			 *
			 * @type {Scene}
			 */

		this.scene = scene;

		/**
			 * Indicates whether the frustum culling is active.
			 *
			 * @type {Boolean}
			 * @default false
			 */

		this.enabled = false;

		/**
			 * A camera.
			 *
			 * @type {PerspectiveCamera}
			 */

		this.cullCamera = new three.PerspectiveCamera(20, 1.77, 0.5, 5);
		this.cullCamera.matrixAutoUpdate = false;

		/**
			 * A spherical coordinate system.
			 *
			 * @type {Spherical}
			 */

		this.s = new three.Spherical(5, Math.PI / 3, Math.PI * 1.75);

		/**
			 * A delta time.
			 *
			 * @type {String}
			 */

		this.delta = "";

		/**
			 * A point cloud that visualises the culled octants.
			 *
			 * @type {Points}
			 */

		this.culledOctants = new three.Points(
			new three.BufferGeometry(),
			new three.PointsMaterial({
				color: 0xccff00,
				sizeAttenuation: false,
				size: 2
			})
		);

		this.culledOctants.visible = false;

		/**
			 * A camera helper.
			 *
			 * @type {CameraHelper}
			 */

		this.cameraHelper = new three.CameraHelper(this.cullCamera);
		this.cameraHelper.visible = false;

	};

	/**
		 * Updates the cull camera.
		 *
		 * @private
		 */

	FrustumCuller.prototype.updateCamera = function updateCamera () {

		var cullCamera = this.cullCamera;

		cullCamera.position.setFromSpherical(this.s);
		cullCamera.lookAt(this.scene.position);

		cullCamera.updateMatrix();
		cullCamera.updateMatrixWorld();

		frustum.setFromMatrix(
			matrix4.multiplyMatrices(
				cullCamera.projectionMatrix,
				cullCamera.matrixWorldInverse
			)
		);

	};

	/**
		 * Culls the octree.
		 */

	FrustumCuller.prototype.cull = function cull () {

		var culledOctants = this.culledOctants;

		var t0;
		var i, j, l;
		var octant, octants;
		var positions;

		if(this.enabled) {

			this.updateCamera();

			t0 = performance.now();
			octants = this.octree.cull(frustum);

			this.delta = (performance.now() - t0).toFixed(2) + " ms";

			if(octants.length > 0) {

				positions = new Float32Array(octants.length * 3 * 2);

				for(i = 0, j = 0, l = octants.length; i < l; ++i) {

					octant = octants[i];
					positions[j++] = octant.min.x;
					positions[j++] = octant.min.y;
					positions[j++] = octant.min.z;
					positions[j++] = octant.max.x;
					positions[j++] = octant.max.y;
					positions[j++] = octant.max.z;

				}

				culledOctants.geometry.removeAttribute("position");
				culledOctants.geometry.addAttribute("position", new three.BufferAttribute(positions, 3));

				this.scene.remove(culledOctants);
				this.scene.add(culledOctants);

			} else {

				this.scene.remove(culledOctants);

			}

			this.cameraHelper.update();

		}

	};

	/**
		 * Registers configuration options.
		 *
		 * @param {GUI} gui - A GUI.
		 */

	FrustumCuller.prototype.configure = function configure (gui) {
			var this$1 = this;


		var folder = gui.addFolder("Frustum Culling");

		folder.add(this, "enabled").onChange(function () {

			this$1.cameraHelper.visible = this$1.culledOctants.visible = this$1.enabled;
			this$1.cull();

		});

		folder.add(this, "delta").listen();
		folder.open();

		var subFolder = folder.addFolder("Camera Adjustment");

		subFolder.add(this.s, "radius").min(0.1).max(10.0).step(0.1).onChange(function () {

			this$1.cull();

		});

		subFolder.add(this.s, "phi").min(1e-6).max(Math.PI - 1e-6).onChange(function () {

			this$1.cull();

		});

		subFolder.add(this.s, "theta").min(0.0).max(Math.PI * 2.0).onChange(function () {

			this$1.cull();

		});

	};

	/**
	 * A point octree demo application.
	 */

	var PointOctreeDemo = (function (Demo$$1) {
		function PointOctreeDemo(renderer) {

			Demo$$1.call(this, renderer);

			/**
			 * A point cloud.
			 *
			 * @type {Points}
			 * @private
			 */

			this.points = null;

			/**
			 * An octree helper.
			 *
			 * @type {OctreeHelper}
			 * @private
			 */

			this.octreeHelper = null;

			/**
			 * An octree raycaster.
			 *
			 * @type {OctreeRaycaster}
			 * @private
			 */

			this.octreeRaycaster = null;

			/**
			 * A frustum culler.
			 *
			 * @type {FrustumCuller}
			 * @private
			 */

			this.frustumCuller = null;

		}

		if ( Demo$$1 ) PointOctreeDemo.__proto__ = Demo$$1;
		PointOctreeDemo.prototype = Object.create( Demo$$1 && Demo$$1.prototype );
		PointOctreeDemo.prototype.constructor = PointOctreeDemo;

		/**
		 * Creates the scene.
		 */

		PointOctreeDemo.prototype.initialise = function initialise () {

			var scene = this.scene;
			var camera = this.camera;
			var renderer = this.renderer;

			// Fog.

			scene.fog = new three.FogExp2(0x0d0d0d, 0.025);

			// Controls.

			this.controls = new three.OrbitControls(camera, renderer.domElement);
			this.controls.maxDistance = 60;

			// Camera.

			camera.near = 0.1;
			camera.far = 200;
			camera.position.set(10, 6, 10);
			camera.lookAt(this.controls.target);

			// Points.

			var points = (function generatePoints() {

				function createPlaneGeometry(particles, n, zBase, zBias) {

					var geometry = new three.BufferGeometry();
					var positions = new Float32Array(particles * 3);
					var n2 = n / 2;

					var x, y, z, i, l;

					for(i = 0, l = positions.length; i < l; i += 3) {

						x = Math.random() * n - n2;
						y = Math.random() * n - n2;
						z = zBase + (Math.random() * zBias * 2 - zBias);

						positions[i] = x;
						positions[i + 1] = y;
						positions[i + 2] = z;

					}

					geometry.addAttribute("position", new three.BufferAttribute(positions, 3));

					return geometry;

				}

				var points = new three.Object3D();

				var w = 256;
				var h = 256;

				var d = 8;

				var size = 6;
				var zStep = size / (d - 1);

				var z = size * -0.5;
				var p;

				var material = new three.PointsMaterial({
					color: 0xc00000,
					sizeAttenuation: false,
					size: 1
				});

				console.log("Generating", w * h * d, "points...");

				while(d-- > 0) {

					p = new three.Points(createPlaneGeometry(w * h, size, z, 0.25), material);
					material = material.clone();
					z += zStep;

					points.add(p);

				}

				return points;

			}());

			this.points = points;
			scene.add(points);

			// Octree.

			var octree = (function createOctree(points) {

				var v = new three.Vector3();
				var bbox = new three.Box3();
				bbox.setFromObject(scene);

				var t0 = performance.now();

				var d, p, i, l;
				var array;

				var octree = new PointOctree(bbox.min, bbox.max, 0.0, 8, 5);

				for(d = points.children.length - 1; d >= 0; --d) {

					p = points.children[d];
					array = p.geometry.getAttribute("position").array;

					for(i = 0, l = array.length; i < l; i += 3) {

						octree.put(v.fromArray(array, i), p);

					}

				}

				console.log("Octree:", octree, "created in", (performance.now() - t0).toFixed(2) + " ms");

				return octree;

			}(points));

			// Octree Helper.

			var octreeHelper = (function createOctreeHelper(octree) {

				var t0 = performance.now();
				var octreeHelper = new OctreeHelper(octree);
				octreeHelper.visible = false;

				console.log("OctreeHelper:", octreeHelper, "created in", (performance.now() - t0).toFixed(2) + " ms");

				return octreeHelper;

			}(octree));

			this.octreeHelper = octreeHelper;
			scene.add(octreeHelper);

			// Raycasting.

			this.raycaster = new OctreeRaycaster$1(octree, camera, points);

			renderer.domElement.parentNode.addEventListener("mousemove", this.raycaster);
			scene.add(this.raycaster.selectedPoint);

			// Frustum culling.

			this.frustumCuller = new FrustumCuller(octree, scene);

			scene.add(this.frustumCuller.cameraHelper);

		};

		/**
		 * Renders this demo.
		 *
		 * @param {Number} delta - The time since the last frame in seconds.
		 */

		PointOctreeDemo.prototype.render = function render (delta) {

			this.renderer.render(this.scene, this.camera);

		};

		/**
		 * Registers configuration options.
		 *
		 * @param {GUI} gui - A GUI.
		 */

		PointOctreeDemo.prototype.configure = function configure (gui) {

			var points = this.points;
			var octreeHelper = this.octreeHelper;

			this.raycaster.configure(gui);
			this.frustumCuller.configure(gui);

			var params = {
				"level mask": octreeHelper.children.length
			};

			var folder = gui.addFolder("Points");
			folder.add(points, "visible");
			folder.open();

			folder = gui.addFolder("Octree Helper");
			folder.add(octreeHelper, "visible");

			folder.add(params, "level mask").min(0).max(octreeHelper.children.length).step(1).onChange(function() {

				var i, l;

				for(i = 0, l = octreeHelper.children.length; i < l; ++i) {

					octreeHelper.children[i].visible = (params["level mask"] === octreeHelper.children.length || i === params["level mask"]);

				}

			});

			folder.open();

		};

		return PointOctreeDemo;
	}(Demo));

	/**
	 * A demo application.
	 */

	var App = function App() {

		/**
			 * A clock.
			 *
			 * @type {Clock}
			 * @private
			 */

		this.clock = new three.Clock();

		/**
			 * A renderer.
			 *
			 * @type {WebGLRenderer}
			 * @private
			 */

		this.renderer = new three.WebGLRenderer({
			logarithmicDepthBuffer: true,
			antialias: true
		});

		this.renderer.setSize(window.innerWidth, window.innerHeight);
		this.renderer.setClearColor(0x000000);
		this.renderer.setPixelRatio(window.devicePixelRatio);

		/**
			 * Statistics.
			 *
			 * @type {Stats}
			 * @private
			 */

		this.stats = (function() {

			var stats = new Stats();
			stats.showPanel(0);
			stats.dom.id = "stats";

			return stats;

		}());

		/**
			 * Available demos.
			 *
			 * @type {Map}
			 * @private
			 */

		this.demos = (function(renderer) {

			var demos = new Map();

			demos.set("point-octree", new PointOctreeDemo(renderer));

			return demos;

		}(this.renderer));

		/**
			 * The key of the current demo.
			 *
			 * @type {String}
			 * @private
			 */

		this.key = (function(demos) {

			var key = window.location.hash.slice(1);

			if(key.length === 0 || !demos.has(key)) {

				key = demos.keys().next().value;

			}

			return key;

		}(this.demos));

	};

	/**
		 * Initialises the demo.
		 *
		 * @param {HTMLElement} viewport - The viewport.
		 * @param {HTMLElement} aside - A secondary DOM container.
		 * @param {HTMLElement} loadingMessage - A loading message.
		 */

	App.prototype.initialise = function initialise (viewport, aside, loadingMessage) {

		var app = this;

		var renderer = this.renderer;
		var clock = this.clock;
		var stats = this.stats;
		var demos = this.demos;

		var demo = null;
		var gui = null;

		viewport.appendChild(renderer.domElement);
		aside.appendChild(stats.dom);

		/**
			 * Activates the currently selected demo.
			 *
			 * @private
			 */

		function activateDemo() {

			demo.initialise();

			demo.camera.aspect = window.innerWidth / window.innerHeight;
			demo.camera.updateProjectionMatrix();

			gui = new dat.GUI({ autoPlace: false });
			gui.add(app, "key", Array.from(demos.keys())).onChange(loadDemo);
			demo.configure(gui);
			aside.appendChild(gui.domElement);

			loadingMessage.style.display = "none";
			renderer.domElement.style.visibility = "visible";

		}

		/**
			 * Loads the currently selected demo.
			 *
			 * @private
			 */

		function loadDemo() {

			var size = renderer.getSize();

			loadingMessage.style.display = "block";
			renderer.domElement.style.visibility = "hidden";

			if(gui !== null) {

				gui.destroy();
				aside.removeChild(gui.domElement);

			}

			if(demo !== null) {

				demo.reset();
				renderer.setSize(size.width, size.height);

			}

			demo = demos.get(app.key);
			demo.load(activateDemo);

			// Update the url.
			window.location.hash = app.key;

		}

		loadDemo();

		/**
			 * Toggles the visibility of the interface on alt key press.
			 *
			 * @private
			 * @param {Event} event - An event.
			 */

		document.addEventListener("keydown", function onKeyDown(event) {

			if(event.altKey) {

				event.preventDefault();
				aside.style.visibility = (aside.style.visibility === "hidden") ? "visible" : "hidden";

			}

		});

		/**
			 * Handles browser resizing.
			 *
			 * @private
			 * @param {Event} event - An event.
			 */

		window.addEventListener("resize", (function() {

			var id = 0;

			function handleResize(event) {

				var width = event.target.innerWidth;
				var height = event.target.innerHeight;

				renderer.setSize(width, height);
				demo.camera.aspect = width / height;
				demo.camera.updateProjectionMatrix();

				id = 0;

			}

			return function onResize(event) {

				if(id === 0) {

					id = setTimeout(handleResize, 66, event);

				}

			};

		}()));

		/**
			 * The main render loop.
			 *
			 * @private
			 * @param {DOMHighResTimeStamp} now - An execution timestamp.
			 */

		(function render(now) {

			var delta = clock.getDelta();

			requestAnimationFrame(render);

			stats.begin();

			demo.render(delta);

			stats.end();

		}());

	};

	/**
	 * Starts the program.
	 *
	 * @private
	 * @param {Event} event - An event.
	 */

	window.addEventListener("load", function main(event) {

		var viewport = document.getElementById("viewport");
		var loadingMessage = viewport.children[0];
		var aside = document.getElementById("aside");

		var app = new App();

		window.removeEventListener("load", main);
		aside.style.visibility = "visible";

		app.initialise(viewport, aside, loadingMessage);

	});

}(THREE,dat,Stats));
