/**
 * sparse-octree v5.0.1 build Thu Jul 05 2018
 * https://github.com/vanruesc/sparse-octree
 * Copyright 2018 Raoul van Rüschen, Zlib
 */
(function (global, factory) {
  typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports, require('three')) :
  typeof define === 'function' && define.amd ? define(['exports', 'three'], factory) :
  (factory((global.SPARSEOCTREE = {}),global.THREE));
}(this, (function (exports,three) { 'use strict';

  var classCallCheck = function (instance, Constructor) {
    if (!(instance instanceof Constructor)) {
      throw new TypeError("Cannot call a class as a function");
    }
  };

  var createClass = function () {
    function defineProperties(target, props) {
      for (var i = 0; i < props.length; i++) {
        var descriptor = props[i];
        descriptor.enumerable = descriptor.enumerable || false;
        descriptor.configurable = true;
        if ("value" in descriptor) descriptor.writable = true;
        Object.defineProperty(target, descriptor.key, descriptor);
      }
    }

    return function (Constructor, protoProps, staticProps) {
      if (protoProps) defineProperties(Constructor.prototype, protoProps);
      if (staticProps) defineProperties(Constructor, staticProps);
      return Constructor;
    };
  }();

  var get = function get(object, property, receiver) {
    if (object === null) object = Function.prototype;
    var desc = Object.getOwnPropertyDescriptor(object, property);

    if (desc === undefined) {
      var parent = Object.getPrototypeOf(object);

      if (parent === null) {
        return undefined;
      } else {
        return get(parent, property, receiver);
      }
    } else if ("value" in desc) {
      return desc.value;
    } else {
      var getter = desc.get;

      if (getter === undefined) {
        return undefined;
      }

      return getter.call(receiver);
    }
  };

  var inherits = function (subClass, superClass) {
    if (typeof superClass !== "function" && superClass !== null) {
      throw new TypeError("Super expression must either be null or a function, not " + typeof superClass);
    }

    subClass.prototype = Object.create(superClass && superClass.prototype, {
      constructor: {
        value: subClass,
        enumerable: false,
        writable: true,
        configurable: true
      }
    });
    if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass;
  };

  var possibleConstructorReturn = function (self, call) {
    if (!self) {
      throw new ReferenceError("this hasn't been initialised - super() hasn't been called");
    }

    return call && (typeof call === "object" || typeof call === "function") ? call : self;
  };

  var toConsumableArray = function (arr) {
    if (Array.isArray(arr)) {
      for (var i = 0, arr2 = Array(arr.length); i < arr.length; i++) arr2[i] = arr[i];

      return arr2;
    } else {
      return Array.from(arr);
    }
  };

  var c = new three.Vector3();

  var Octant = function () {
  	function Octant() {
  		var min = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : new three.Vector3();
  		var max = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : new three.Vector3();
  		classCallCheck(this, Octant);


  		this.min = min;

  		this.max = max;

  		this.children = null;
  	}

  	createClass(Octant, [{
  		key: "getCenter",
  		value: function getCenter() {
  			var target = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : new three.Vector3();


  			return target.addVectors(this.min, this.max).multiplyScalar(0.5);
  		}
  	}, {
  		key: "getDimensions",
  		value: function getDimensions() {
  			var target = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : new three.Vector3();


  			return target.subVectors(this.max, this.min);
  		}
  	}, {
  		key: "split",
  		value: function split() {

  			var min = this.min;
  			var max = this.max;
  			var mid = this.getCenter(c);

  			var children = this.children = [null, null, null, null, null, null, null, null];

  			var i = void 0,
  			    combination = void 0;

  			for (i = 0; i < 8; ++i) {

  				combination = pattern[i];

  				children[i] = new this.constructor(new three.Vector3(combination[0] === 0 ? min.x : mid.x, combination[1] === 0 ? min.y : mid.y, combination[2] === 0 ? min.z : mid.z), new three.Vector3(combination[0] === 0 ? mid.x : max.x, combination[1] === 0 ? mid.y : max.y, combination[2] === 0 ? mid.z : max.z));
  			}
  		}
  	}]);
  	return Octant;
  }();

  var pattern = [new Uint8Array([0, 0, 0]), new Uint8Array([0, 0, 1]), new Uint8Array([0, 1, 0]), new Uint8Array([0, 1, 1]), new Uint8Array([1, 0, 0]), new Uint8Array([1, 0, 1]), new Uint8Array([1, 1, 0]), new Uint8Array([1, 1, 1])];

  var edges = [new Uint8Array([0, 4]), new Uint8Array([1, 5]), new Uint8Array([2, 6]), new Uint8Array([3, 7]), new Uint8Array([0, 2]), new Uint8Array([1, 3]), new Uint8Array([4, 6]), new Uint8Array([5, 7]), new Uint8Array([0, 1]), new Uint8Array([2, 3]), new Uint8Array([4, 5]), new Uint8Array([6, 7])];

  var c$1 = new three.Vector3();

  var CubicOctant = function () {
  	function CubicOctant() {
  		var min = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : new three.Vector3();
  		var size = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 0;
  		classCallCheck(this, CubicOctant);


  		this.min = min;

  		this.size = size;

  		this.children = null;
  	}

  	createClass(CubicOctant, [{
  		key: "getCenter",
  		value: function getCenter() {
  			var target = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : new three.Vector3();


  			return target.copy(this.min).addScalar(this.size * 0.5);
  		}
  	}, {
  		key: "getDimensions",
  		value: function getDimensions() {
  			var target = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : new three.Vector3();


  			return target.set(this.size, this.size, this.size);
  		}
  	}, {
  		key: "split",
  		value: function split() {

  			var min = this.min;
  			var mid = this.getCenter(c$1);
  			var halfSize = this.size * 0.5;

  			var children = this.children = [null, null, null, null, null, null, null, null];

  			var i = void 0,
  			    combination = void 0;

  			for (i = 0; i < 8; ++i) {

  				combination = pattern[i];

  				children[i] = new this.constructor(new three.Vector3(combination[0] === 0 ? min.x : mid.x, combination[1] === 0 ? min.y : mid.y, combination[2] === 0 ? min.z : mid.z), halfSize);
  			}
  		}
  	}, {
  		key: "max",
  		get: function get$$1() {

  			return this.min.clone().addScalar(this.size);
  		}
  	}]);
  	return CubicOctant;
  }();

  var IteratorResult = function () {
  	function IteratorResult() {
  		var value = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : null;
  		var done = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : false;
  		classCallCheck(this, IteratorResult);


  		this.value = value;

  		this.done = done;
  	}

  	createClass(IteratorResult, [{
  		key: "reset",
  		value: function reset() {

  			this.value = null;
  			this.done = false;
  		}
  	}]);
  	return IteratorResult;
  }();

  var b = new three.Box3();

  var OctantIterator = function () {
  	function OctantIterator(octree) {
  		var region = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : null;
  		classCallCheck(this, OctantIterator);


  		this.octree = octree;

  		this.region = region;

  		this.cull = region !== null;

  		this.result = new IteratorResult();

  		this.trace = null;

  		this.indices = null;

  		this.reset();
  	}

  	createClass(OctantIterator, [{
  		key: "reset",
  		value: function reset() {

  			var root = this.octree.root;

  			this.trace = [];
  			this.indices = [];

  			if (root !== null) {

  				b.min = root.min;
  				b.max = root.max;

  				if (!this.cull || this.region.intersectsBox(b)) {

  					this.trace.push(root);
  					this.indices.push(0);
  				}
  			}

  			this.result.reset();

  			return this;
  		}
  	}, {
  		key: "next",
  		value: function next() {

  			var cull = this.cull;
  			var region = this.region;
  			var indices = this.indices;
  			var trace = this.trace;

  			var octant = null;
  			var depth = trace.length - 1;

  			var index = void 0,
  			    children = void 0,
  			    child = void 0;

  			while (octant === null && depth >= 0) {

  				index = indices[depth]++;
  				children = trace[depth].children;

  				if (index < 8) {

  					if (children !== null) {

  						child = children[index];

  						if (cull) {

  							b.min = child.min;
  							b.max = child.max;

  							if (!region.intersectsBox(b)) {
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
  			this.result.done = octant === null;

  			return this.result;
  		}
  	}, {
  		key: "return",
  		value: function _return(value) {

  			this.result.value = value;
  			this.result.done = true;

  			return this.result;
  		}
  	}, {
  		key: Symbol.iterator,
  		value: function value() {

  			return this;
  		}
  	}]);
  	return OctantIterator;
  }();

  var v = [new three.Vector3(), new three.Vector3(), new three.Vector3()];

  var b$1 = new three.Box3();

  var r = new three.Ray();

  var octantTable = [new Uint8Array([4, 2, 1]), new Uint8Array([5, 3, 8]), new Uint8Array([6, 8, 3]), new Uint8Array([7, 8, 8]), new Uint8Array([8, 6, 5]), new Uint8Array([8, 7, 8]), new Uint8Array([8, 8, 7]), new Uint8Array([8, 8, 8])];

  var flags = 0;

  function findEntryOctant(tx0, ty0, tz0, txm, tym, tzm) {

  	var entry = 0;

  	if (tx0 > ty0 && tx0 > tz0) {
  		if (tym < tx0) {

  			entry |= 2;
  		}

  		if (tzm < tx0) {

  			entry |= 1;
  		}
  	} else if (ty0 > tz0) {
  		if (txm < ty0) {

  			entry |= 4;
  		}

  		if (tzm < ty0) {

  			entry |= 1;
  		}
  	} else {
  		if (txm < tz0) {

  			entry |= 4;
  		}

  		if (tym < tz0) {

  			entry |= 2;
  		}
  	}

  	return entry;
  }

  function findNextOctant(currentOctant, tx1, ty1, tz1) {

  	var min = void 0;
  	var exit = 0;

  	if (tx1 < ty1) {

  		min = tx1;
  		exit = 0;
  	} else {

  		min = ty1;
  		exit = 1;
  	}

  	if (tz1 < min) {

  		exit = 2;
  	}

  	return octantTable[currentOctant][exit];
  }

  function raycastOctant(octant, tx0, ty0, tz0, tx1, ty1, tz1, raycaster, intersects) {

  	var children = octant.children;

  	var currentOctant = void 0;
  	var txm = void 0,
  	    tym = void 0,
  	    tzm = void 0;

  	if (tx1 >= 0.0 && ty1 >= 0.0 && tz1 >= 0.0) {

  		if (children === null) {
  			intersects.push(octant);
  		} else {
  			txm = 0.5 * (tx0 + tx1);
  			tym = 0.5 * (ty0 + ty1);
  			tzm = 0.5 * (tz0 + tz1);

  			currentOctant = findEntryOctant(tx0, ty0, tz0, txm, tym, tzm);

  			do {

  				switch (currentOctant) {

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

  						currentOctant = 8;
  						break;

  				}
  			} while (currentOctant < 8);
  		}
  	}
  }

  var OctreeRaycaster = function () {
  	function OctreeRaycaster() {
  		classCallCheck(this, OctreeRaycaster);
  	}

  	createClass(OctreeRaycaster, null, [{
  		key: "intersectOctree",
  		value: function intersectOctree(octree, raycaster, intersects) {
  			var min = b$1.min.set(0, 0, 0);
  			var max = b$1.max.subVectors(octree.max, octree.min);

  			var dimensions = octree.getDimensions(v[0]);
  			var halfDimensions = v[1].copy(dimensions).multiplyScalar(0.5);

  			var origin = r.origin.copy(raycaster.ray.origin);
  			var direction = r.direction.copy(raycaster.ray.direction);

  			var invDirX = void 0,
  			    invDirY = void 0,
  			    invDirZ = void 0;
  			var tx0 = void 0,
  			    tx1 = void 0,
  			    ty0 = void 0,
  			    ty1 = void 0,
  			    tz0 = void 0,
  			    tz1 = void 0;

  			origin.sub(octree.getCenter(v[2])).add(halfDimensions);

  			flags = 0;

  			if (direction.x < 0.0) {

  				origin.x = dimensions.x - origin.x;
  				direction.x = -direction.x;
  				flags |= 4;
  			}

  			if (direction.y < 0.0) {

  				origin.y = dimensions.y - origin.y;
  				direction.y = -direction.y;
  				flags |= 2;
  			}

  			if (direction.z < 0.0) {

  				origin.z = dimensions.z - origin.z;
  				direction.z = -direction.z;
  				flags |= 1;
  			}

  			invDirX = 1.0 / direction.x;
  			invDirY = 1.0 / direction.y;
  			invDirZ = 1.0 / direction.z;

  			tx0 = (min.x - origin.x) * invDirX;
  			tx1 = (max.x - origin.x) * invDirX;
  			ty0 = (min.y - origin.y) * invDirY;
  			ty1 = (max.y - origin.y) * invDirY;
  			tz0 = (min.z - origin.z) * invDirZ;
  			tz1 = (max.z - origin.z) * invDirZ;

  			if (Math.max(Math.max(tx0, ty0), tz0) < Math.min(Math.min(tx1, ty1), tz1)) {
  				raycastOctant(octree.root, tx0, ty0, tz0, tx1, ty1, tz1, raycaster, intersects);
  			}
  		}
  	}]);
  	return OctreeRaycaster;
  }();

  var b$2 = new three.Box3();

  function _getDepth(octant) {

  	var children = octant.children;

  	var result = 0;
  	var i = void 0,
  	    l = void 0,
  	    d = void 0;

  	if (children !== null) {

  		for (i = 0, l = children.length; i < l; ++i) {

  			d = 1 + _getDepth(children[i]);

  			if (d > result) {

  				result = d;
  			}
  		}
  	}

  	return result;
  }

  function _cull(octant, region, result) {

  	var children = octant.children;

  	var i = void 0,
  	    l = void 0;

  	b$2.min = octant.min;
  	b$2.max = octant.max;

  	if (region.intersectsBox(b$2)) {

  		if (children !== null) {

  			for (i = 0, l = children.length; i < l; ++i) {

  				_cull(children[i], region, result);
  			}
  		} else {

  			result.push(octant);
  		}
  	}
  }

  function _findOctantsByLevel(octant, level, depth, result) {

  	var children = octant.children;

  	var i = void 0,
  	    l = void 0;

  	if (depth === level) {

  		result.push(octant);
  	} else if (children !== null) {

  		++depth;

  		for (i = 0, l = children.length; i < l; ++i) {

  			_findOctantsByLevel(children[i], level, depth, result);
  		}
  	}
  }

  var Octree = function () {
  	function Octree(min, max) {
  		classCallCheck(this, Octree);


  		this.root = min !== undefined && max !== undefined ? new Octant(min, max) : null;
  	}

  	createClass(Octree, [{
  		key: "getCenter",
  		value: function getCenter(target) {

  			return this.root.getCenter(target);
  		}
  	}, {
  		key: "getDimensions",
  		value: function getDimensions(target) {

  			return this.root.getDimensions(target);
  		}
  	}, {
  		key: "getDepth",
  		value: function getDepth() {

  			return _getDepth(this.root);
  		}
  	}, {
  		key: "cull",
  		value: function cull(region) {

  			var result = [];

  			_cull(this.root, region, result);

  			return result;
  		}
  	}, {
  		key: "findOctantsByLevel",
  		value: function findOctantsByLevel(level) {

  			var result = [];

  			_findOctantsByLevel(this.root, level, 0, result);

  			return result;
  		}
  	}, {
  		key: "raycast",
  		value: function raycast(raycaster) {
  			var intersects = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : [];


  			OctreeRaycaster.intersectOctree(this, raycaster, intersects);

  			return intersects;
  		}
  	}, {
  		key: "leaves",
  		value: function leaves(region) {

  			return new OctantIterator(this, region);
  		}
  	}, {
  		key: Symbol.iterator,
  		value: function value() {

  			return new OctantIterator(this);
  		}
  	}, {
  		key: "min",
  		get: function get$$1() {

  			return this.root.min;
  		}
  	}, {
  		key: "max",
  		get: function get$$1() {

  			return this.root.max;
  		}
  	}, {
  		key: "children",
  		get: function get$$1() {

  			return this.root.children;
  		}
  	}]);
  	return Octree;
  }();

  var p = new three.Vector3();

  var PointOctant = function (_Octant) {
  	inherits(PointOctant, _Octant);

  	function PointOctant(min, max) {
  		classCallCheck(this, PointOctant);

  		var _this = possibleConstructorReturn(this, (PointOctant.__proto__ || Object.getPrototypeOf(PointOctant)).call(this, min, max));

  		_this.points = null;

  		_this.data = null;

  		return _this;
  	}

  	createClass(PointOctant, [{
  		key: "distanceToSquared",
  		value: function distanceToSquared(point) {

  			var clampedPoint = p.copy(point).clamp(this.min, this.max);

  			return clampedPoint.sub(point).lengthSquared();
  		}
  	}, {
  		key: "distanceToCenterSquared",
  		value: function distanceToCenterSquared(point) {

  			var center = this.getCenter(p);

  			var dx = point.x - center.x;
  			var dy = point.y - center.x;
  			var dz = point.z - center.z;

  			return dx * dx + dy * dy + dz * dz;
  		}
  	}, {
  		key: "contains",
  		value: function contains(point, bias) {

  			var min = this.min;
  			var max = this.max;

  			return point.x >= min.x - bias && point.y >= min.y - bias && point.z >= min.z - bias && point.x <= max.x + bias && point.y <= max.y + bias && point.z <= max.z + bias;
  		}
  	}, {
  		key: "redistribute",
  		value: function redistribute(bias) {

  			var children = this.children;
  			var points = this.points;
  			var data = this.data;

  			var i = void 0,
  			    j = void 0,
  			    il = void 0,
  			    jl = void 0;
  			var child = void 0,
  			    point = void 0,
  			    entry = void 0;

  			if (children !== null && points !== null) {

  				for (i = 0, il = points.length; i < il; ++i) {

  					point = points[i];
  					entry = data[i];

  					for (j = 0, jl = children.length; j < jl; ++j) {

  						child = children[j];

  						if (child.contains(point, bias)) {

  							if (child.points === null) {

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
  		}
  	}, {
  		key: "merge",
  		value: function merge() {

  			var children = this.children;

  			var i = void 0,
  			    l = void 0;
  			var child = void 0;

  			if (children !== null) {

  				this.points = [];
  				this.data = [];

  				for (i = 0, l = children.length; i < l; ++i) {

  					child = children[i];

  					if (child.points !== null) {
  						var _points, _data;

  						(_points = this.points).push.apply(_points, toConsumableArray(child.points));
  						(_data = this.data).push.apply(_data, toConsumableArray(child.data));
  					}
  				}

  				this.children = null;
  			}
  		}
  	}]);
  	return PointOctant;
  }(Octant);

  var RayPointIntersection = function RayPointIntersection(distance, distanceToRay, point) {
  	var object = arguments.length > 3 && arguments[3] !== undefined ? arguments[3] : null;
  	classCallCheck(this, RayPointIntersection);


  	this.distance = distance;

  	this.distanceToRay = distanceToRay;

  	this.point = point;

  	this.object = object;
  };

  var THRESHOLD = 1e-6;

  function _countPoints(octant) {

  	var children = octant.children;

  	var result = 0;
  	var i = void 0,
  	    l = void 0;

  	if (children !== null) {

  		for (i = 0, l = children.length; i < l; ++i) {

  			result += _countPoints(children[i]);
  		}
  	} else if (octant.points !== null) {

  		result = octant.points.length;
  	}

  	return result;
  }

  function _put(point, data, octree, octant, depth) {

  	var children = octant.children;
  	var exists = false;
  	var done = false;
  	var i = void 0,
  	    l = void 0;

  	if (octant.contains(point, octree.bias)) {

  		if (children === null) {

  			if (octant.points === null) {

  				octant.points = [];
  				octant.data = [];
  			} else {

  				for (i = 0, l = octant.points.length; !exists && i < l; ++i) {

  					exists = octant.points[i].equals(point);
  				}
  			}

  			if (exists) {

  				octant.data[i - 1] = data;
  				done = true;
  			} else if (octant.points.length < octree.maxPoints || depth === octree.maxDepth) {

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

  		if (children !== null) {

  			++depth;

  			for (i = 0, l = children.length; !done && i < l; ++i) {

  				done = _put(point, data, octree, children[i], depth);
  			}
  		}
  	}

  	return done;
  }

  function _remove(point, octree, octant, parent) {

  	var children = octant.children;

  	var result = null;

  	var i = void 0,
  	    l = void 0;
  	var points = void 0,
  	    data = void 0,
  	    last = void 0;

  	if (octant.contains(point, octree.bias)) {

  		if (children !== null) {

  			for (i = 0, l = children.length; result === null && i < l; ++i) {

  				result = _remove(point, octree, children[i], octant);
  			}
  		} else if (octant.points !== null) {

  			points = octant.points;
  			data = octant.data;

  			for (i = 0, l = points.length; i < l; ++i) {

  				if (points[i].equals(point)) {

  					last = l - 1;
  					result = data[i];

  					if (i < last) {
  						points[i] = points[last];
  						data[i] = data[last];
  					}

  					points.pop();
  					data.pop();

  					--octree.pointCount;

  					if (parent !== null && _countPoints(parent) <= octree.maxPoints) {

  						parent.merge();
  					}

  					break;
  				}
  			}
  		}
  	}

  	return result;
  }

  function _fetch(point, octree, octant) {

  	var children = octant.children;

  	var result = null;

  	var i = void 0,
  	    l = void 0;
  	var points = void 0;

  	if (octant.contains(point, octree.bias)) {

  		if (children !== null) {

  			for (i = 0, l = children.length; result === null && i < l; ++i) {

  				result = _fetch(point, octree, children[i]);
  			}
  		} else if (octant.points !== null) {

  			points = octant.points;

  			for (i = 0, l = points.length; result === null && i < l; ++i) {

  				if (point.distanceToSquared(points[i]) <= THRESHOLD) {

  					result = octant.data[i];
  				}
  			}
  		}
  	}

  	return result;
  }

  function _move(point, position, octree, octant, parent, depth) {

  	var children = octant.children;

  	var result = null;

  	var i = void 0,
  	    l = void 0;
  	var points = void 0;

  	if (octant.contains(point, octree.bias)) {

  		if (octant.contains(position, octree.bias)) {
  			if (children !== null) {

  				++depth;

  				for (i = 0, l = children.length; result === null && i < l; ++i) {

  					result = _move(point, position, octree, children[i], octant, depth);
  				}
  			} else if (octant.points !== null) {
  				points = octant.points;

  				for (i = 0, l = points.length; i < l; ++i) {

  					if (point.distanceToSquared(points[i]) <= THRESHOLD) {
  						points[i].copy(position);
  						result = octant.data[i];

  						break;
  					}
  				}
  			}
  		} else {
  			result = _remove(point, octree, octant, parent);

  			_put(position, result, octree, parent, depth - 1);
  		}
  	}

  	return result;
  }

  function _findNearestPoint(point, maxDistance, skipSelf, octant) {

  	var points = octant.points;
  	var children = octant.children;

  	var result = null;
  	var bestDist = maxDistance;

  	var i = void 0,
  	    l = void 0;
  	var p = void 0,
  	    distSq = void 0;

  	var sortedChildren = void 0;
  	var child = void 0,
  	    childResult = void 0;

  	if (children !== null) {
  		sortedChildren = children.map(function (child) {
  			return {
  				octant: child,
  				distance: child.distanceToCenterSquared(point)
  			};
  		}).sort(function (a, b) {
  			return a.distance - b.distance;
  		});

  		for (i = 0, l = sortedChildren.length; i < l; ++i) {
  			child = sortedChildren[i].octant;

  			if (child.contains(point, bestDist)) {

  				childResult = _findNearestPoint(point, bestDist, skipSelf, child);

  				if (childResult !== null) {

  					distSq = childResult.point.distanceToSquared(point);

  					if ((!skipSelf || distSq > 0.0) && distSq < bestDist) {

  						bestDist = distSq;
  						result = childResult;
  					}
  				}
  			}
  		}
  	} else if (points !== null) {

  		for (i = 0, l = points.length; i < l; ++i) {

  			p = points[i];
  			distSq = point.distanceToSquared(p);

  			if ((!skipSelf || distSq > 0.0) && distSq < bestDist) {

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

  function _findPoints(point, radius, skipSelf, octant, result) {

  	var points = octant.points;
  	var children = octant.children;
  	var rSq = radius * radius;

  	var i = void 0,
  	    l = void 0;

  	var p = void 0,
  	    distSq = void 0;
  	var child = void 0;

  	if (children !== null) {

  		for (i = 0, l = children.length; i < l; ++i) {

  			child = children[i];

  			if (child.contains(point, radius)) {

  				_findPoints(point, radius, skipSelf, child, result);
  			}
  		}
  	} else if (points !== null) {

  		for (i = 0, l = points.length; i < l; ++i) {

  			p = points[i];
  			distSq = point.distanceToSquared(p);

  			if ((!skipSelf || distSq > 0.0) && distSq <= rSq) {

  				result.push({
  					point: p.clone(),
  					data: octant.data[i]
  				});
  			}
  		}
  	}
  }

  var PointOctree = function (_Octree) {
  	inherits(PointOctree, _Octree);

  	function PointOctree(min, max) {
  		var bias = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : 0.0;
  		var maxPoints = arguments.length > 3 && arguments[3] !== undefined ? arguments[3] : 8;
  		var maxDepth = arguments.length > 4 && arguments[4] !== undefined ? arguments[4] : 8;
  		classCallCheck(this, PointOctree);

  		var _this = possibleConstructorReturn(this, (PointOctree.__proto__ || Object.getPrototypeOf(PointOctree)).call(this));

  		_this.root = new PointOctant(min, max);

  		_this.bias = Math.max(0.0, bias);

  		_this.maxPoints = Math.max(1, Math.round(maxPoints));

  		_this.maxDepth = Math.max(0, Math.round(maxDepth));

  		_this.pointCount = 0;

  		return _this;
  	}

  	createClass(PointOctree, [{
  		key: "countPoints",
  		value: function countPoints(octant) {

  			return _countPoints(octant);
  		}
  	}, {
  		key: "put",
  		value: function put(point, data) {

  			return _put(point, data, this, this.root, 0);
  		}
  	}, {
  		key: "remove",
  		value: function remove(point) {

  			return _remove(point, this, this.root, null);
  		}
  	}, {
  		key: "fetch",
  		value: function fetch(point) {

  			return _fetch(point, this, this.root);
  		}
  	}, {
  		key: "move",
  		value: function move(point, position) {

  			return _move(point, position, this, this.root, null, 0);
  		}
  	}, {
  		key: "findNearestPoint",
  		value: function findNearestPoint(point) {
  			var maxDistance = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : Infinity;
  			var skipSelf = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : false;


  			return _findNearestPoint(point, maxDistance, skipSelf, this.root);
  		}
  	}, {
  		key: "findPoints",
  		value: function findPoints(point, radius) {
  			var skipSelf = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : false;


  			var result = [];

  			_findPoints(point, radius, skipSelf, this.root, result);

  			return result;
  		}
  	}, {
  		key: "raycast",
  		value: function raycast(raycaster) {
  			var intersects = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : [];


  			var octants = get(PointOctree.prototype.__proto__ || Object.getPrototypeOf(PointOctree.prototype), "raycast", this).call(this, raycaster);

  			if (octants.length > 0) {
  				this.testPoints(octants, raycaster, intersects);
  			}

  			return intersects;
  		}
  	}, {
  		key: "testPoints",
  		value: function testPoints(octants, raycaster, intersects) {

  			var threshold = raycaster.params.Points.threshold;
  			var thresholdSq = threshold * threshold;

  			var intersectPoint = void 0;
  			var distance = void 0,
  			    distanceToRay = void 0;
  			var rayPointDistanceSq = void 0;

  			var i = void 0,
  			    j = void 0,
  			    il = void 0,
  			    jl = void 0;
  			var octant = void 0,
  			    points = void 0,
  			    point = void 0;

  			for (i = 0, il = octants.length; i < il; ++i) {

  				octant = octants[i];
  				points = octant.points;

  				if (points !== null) {

  					for (j = 0, jl = points.length; j < jl; ++j) {

  						point = points[j];
  						rayPointDistanceSq = raycaster.ray.distanceSqToPoint(point);

  						if (rayPointDistanceSq < thresholdSq) {

  							intersectPoint = raycaster.ray.closestPointToPoint(point, new three.Vector3());
  							distance = raycaster.ray.origin.distanceTo(intersectPoint);

  							if (distance >= raycaster.near && distance <= raycaster.far) {

  								distanceToRay = Math.sqrt(rayPointDistanceSq);

  								intersects.push(new RayPointIntersection(distance, distanceToRay, intersectPoint, octant.data[j]));
  							}
  						}
  					}
  				}
  			}
  		}
  	}]);
  	return PointOctree;
  }(Octree);

  var b$3 = new three.Box3();

  var c$2 = new three.Vector3();

  var u = new three.Vector3();

  var v$1 = new three.Vector3();

  var OctreeUtils = function () {
  	function OctreeUtils() {
  		classCallCheck(this, OctreeUtils);
  	}

  	createClass(OctreeUtils, null, [{
  		key: "recycleOctants",
  		value: function recycleOctants(octant, octants) {

  			var min = octant.min;
  			var mid = octant.getCenter(u);
  			var halfDimensions = octant.getDimensions(v$1).multiplyScalar(0.5);

  			var children = octant.children;
  			var l = octants.length;

  			var i = void 0,
  			    j = void 0;
  			var combination = void 0,
  			    candidate = void 0;

  			for (i = 0; i < 8; ++i) {

  				combination = pattern[i];

  				b$3.min.addVectors(min, c$2.fromArray(combination).multiply(halfDimensions));
  				b$3.max.addVectors(mid, c$2.fromArray(combination).multiply(halfDimensions));

  				for (j = 0; j < l; ++j) {

  					candidate = octants[j];

  					if (candidate !== null && b$3.min.equals(candidate.min) && b$3.max.equals(candidate.max)) {

  						children[i] = candidate;
  						octants[j] = null;

  						break;
  					}
  				}
  			}
  		}
  	}]);
  	return OctreeUtils;
  }();

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
