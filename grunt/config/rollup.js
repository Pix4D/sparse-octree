const resolve = require("rollup-plugin-node-resolve");
const buble = require("rollup-plugin-buble");

module.exports = function(grunt) {

	return {

		options: {
			plugins() {

				return [
					resolve({
						jsnext: true
					})
				].concat(grunt.option("production") ? [buble()] : []);

			}
		},

		lib: {
			options: {
				format: "umd",
				moduleName: "<%= package.name.replace(/-/g, \"\").toUpperCase() %>",
				banner: "<%= banner %>",
				external: ["three"]
			},
			src: "src/index.js",
			dest: "build/<%= package.name %>.js"
		},

		demo: {
			options: {
				globals: {
					"three": "THREE",
					"stats.js": "Stats",
					"dat.gui": "dat"
				},
				external: [
					"three",
					"stats.js",
					"dat.gui"
				],
				format: "iife"
			},
			src: "demo/src/index.js",
			dest: "public/demo/index.js"
		}

	};

};
