/**
 * Isolated 3D Viewer GUI JavaScript
 * Extracted from perry3d.js for standalone use
 * 
 * This file contains the core Three.js setup and interaction code
 * necessary for the 3D molecular structure viewer GUI.
 * 
 * Dependencies:
 * - Three.js library (three.min.js)
 * - TrackballControls.js
 * - CSS2DRenderer.js (for labels)
 */

// ============================================================================
// CONFIGURATION AND SETUP
// ============================================================================

var modelSize = {
    aspect: 2,
    width: 800,
    height: 400
};

var scene = new THREE.Scene();
var camera = new THREE.PerspectiveCamera(30, modelSize.aspect, .5, 200);

// Check for WebGL support
function webglAvailable() {
    try {
        var canvas = document.createElement('canvas');
        return !!(window.WebGLRenderingContext && (
            canvas.getContext('webgl') ||
            canvas.getContext('experimental-webgl')
        ));
    } catch (e) {
        return false;
    }
}

// Initialize renderer based on WebGL availability
var renderer;
if (webglAvailable()) {
    renderer = new THREE.WebGLRenderer({
        preserveDrawingBuffer: true,
        alpha: true
    });
} else {
    renderer = new THREE.CanvasRenderer();
}

// Setup label renderer for 2D text overlays
var labelRenderer = new THREE.CSS2DRenderer();
labelRenderer.domElement.style.position = 'absolute';
labelRenderer.domElement.style.top = '0';
labelRenderer.domElement.style.pointerEvents = 'none';
labelRenderer.domElement.id = 'labelPlace';

// Configure renderer
renderer.setSize(modelSize.width, modelSize.height);
renderer.setClearColor(0xffffff, 0);

// ============================================================================
// CONTROLS SETUP
// ============================================================================

var controls = new THREE.TrackballControls(camera, renderer.domElement);
controls.target.set(0, 0, 0);
controls.position0.set(0, 0, 4); // Initial camera position
controls.rotateSpeed = 3.0;
controls.zoomSpeed = 1.2;
controls.panSpeed = 0.8;
controls.noZoom = false;
controls.noPan = false;
controls.staticMoving = false;
controls.dynamicDampingFactor = 0.15;
controls.reset(); // Place camera at starting position
controls.enabled = false; // Initially disabled until user clicks

// ============================================================================
// LIGHTING SETUP
// ============================================================================

var light = new THREE.DirectionalLight(0xffffff, 1);
scene.add(light);

// Update light position to follow camera
var light_update = function () {
    'use strict';
    light.position.copy(camera.position);
};

light_update();
controls.addEventListener('change', light_update);

// ============================================================================
// RAYCASTING AND MOUSE INTERACTION
// ============================================================================

var mouse = new THREE.Vector2();
var raycaster = new THREE.Raycaster();

// Arrays to store scene objects for interaction
var atoms = [];
var bondsArray = [];

/**
 * Get the first intersected object from mouse click
 * @param {MouseEvent} event - Mouse click event
 * @returns {Object|null} Intersected object or null
 */
function firstIntersectedObject(event) {
    event.preventDefault();
    
    // Calculate mouse position in normalized device coordinates (-1 to +1)
    var container = document.getElementById('model3d_container');
    var rect = container.getBoundingClientRect();
    
    mouse.x = ((event.clientX - rect.left) / modelSize.width) * 2 - 1;
    mouse.y = -((event.clientY - rect.top) / modelSize.height) * 2 + 1;
    
    // Update raycaster with camera and mouse position
    raycaster.setFromCamera(mouse, camera);
    
    // Check intersections with atoms and bonds
    var intersects = raycaster.intersectObjects(atoms);
    var bondIntersects = raycaster.intersectObjects(bondsArray);
    
    // Prioritize atoms over bonds if both are intersected
    if (intersects.length > 0) {
        if (bondIntersects.length > 0) {
            var total = bondIntersects.concat(intersects);
            var tot = [];
            for (var ob of total) {
                tot.push(ob.object);
            }
            var total = raycaster.intersectObjects(tot);
            return total[0];
        }
        return intersects[0];
    }
    
    if (bondIntersects.length > 0) {
        return bondIntersects[0];
    }
    
    return null;
}

// ============================================================================
// MOUSE CLICK HANDLING
// ============================================================================

/**
 * Variable that stores the action to take when an object is clicked
 * This allows for different interaction modes (selection, deletion, etc.)
 */
var takenAction = function () {};

/**
 * Set the active function to execute on mouse click
 * @param {Function} func - Function to execute when object is clicked
 */
function activeFunction(func) {
    takenAction = function (param) {
        func(param);
    };
}

/**
 * Default mouse click handler
 * @param {MouseEvent} event - Mouse click event
 */
function onMouseClick(event) {
    var intersected = firstIntersectedObject(event);
    if (intersected != null) {
        takenAction(intersected.object);
    }
}

// Attach click handler to container
var moleculeClick = onMouseClick;

/**
 * Change the click function handler
 * @param {Function} func - New click handler function
 */
function clickFunction(func) {
    var container = document.getElementById('model3d_container');
    container.removeEventListener('mousedown', moleculeClick, true);
    moleculeClick = func;
    container.addEventListener('mousedown', moleculeClick, true);
}

// ============================================================================
// RENDERING LOOP
// ============================================================================

/**
 * Main rendering loop
 * Updates controls and renders the scene
 */
var render = function () {
    'use strict';
    
    requestAnimationFrame(render);
    controls.update();
    renderer.render(scene, camera);
    labelRenderer.render(scene, camera);
};

// Start the render loop
render();

// ============================================================================
// SCENE MANAGEMENT
// ============================================================================

/**
 * Clear all objects from the scene
 */
function clearScene() {
    for (var i = 0; i < atoms.length; i++) {
        scene.remove(atoms[i]);
    }
    for (var j = 0; j < bondsArray.length; j++) {
        scene.remove(bondsArray[j]);
    }
    
    // Clear any arrows
    while (scene.getObjectByName("arrow")) {
        var arrow = scene.getObjectByName('arrow');
        scene.remove(arrow);
    }
    
    // Clear any orbitals
    while (scene.getObjectByName("orbital")) {
        var orb = scene.getObjectByName('orbital');
        scene.remove(orb);
    }
    
    // Clear any cube files
    while (scene.getObjectByName("cube")) {
        var cube = scene.getObjectByName('cube');
        scene.remove(cube);
    }
    
    bondsArray = [];
    atoms = [];
}

/**
 * Clear selection highlighting
 */
function clearSelection() {
    // This would typically reset selected atom colors
    // Implementation depends on your selection system
}

// ============================================================================
// CANVAS RESIZE HANDLING
// ============================================================================

/**
 * Resize the canvas when container size changes
 * @param {number} new_width - New width for the canvas
 */
var resize_canvas = function(new_width) {
    if (new_width != modelSize.width) {
        modelSize.width = new_width;
        modelSize.height = modelSize.width / modelSize.aspect;
        renderer.setSize(modelSize.width, modelSize.height);
        labelRenderer.setSize(modelSize.width, modelSize.height);
        camera.aspect = modelSize.aspect;
        camera.updateProjectionMatrix();
    }
};

// ============================================================================
// INITIALIZATION
// ============================================================================

/**
 * Initialize the 3D viewer
 * Call this after the DOM is ready and container element exists
 */
function initializeViewer() {
    // Get container element
    var container = document.getElementById('model3d');
    if (!container) {
        console.error('Container element with id "model3d" not found');
        return;
    }
    
    // Append renderer canvas to container
    container.appendChild(renderer.domElement);
    
    // Append label renderer to container
    container.appendChild(labelRenderer.domElement);
    
    // Set canvas ID
    var canvas = renderer.domElement;
    canvas.setAttribute("id", "mol3dContext");
    canvas.style.backgroundColor = "#000000";
    
    // Attach click handler
    var containerElement = document.getElementById('model3d_container');
    if (containerElement) {
        containerElement.addEventListener('mousedown', moleculeClick, true);
    }
    
    // Enable controls when canvas is clicked
    if (container) {
        container.addEventListener('click', function() {
            controls.enabled = true;
        });
    }
    
    // Handle window resize
    window.addEventListener('resize', function() {
        var container = document.getElementById('model3d');
        if (container) {
            resize_canvas(parseInt(container.offsetWidth));
        }
    });
    
    // Initial resize
    resize_canvas(parseInt(container.offsetWidth));
    
    console.log('3D Viewer initialized');
}

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

/**
 * Reset camera to initial position
 */
function resetView() {
    controls.reset();
}

/**
 * Center the camera on a specific point
 * @param {THREE.Vector3} point - Point to center on (optional, defaults to origin)
 */
function centerOnPoint(point) {
    point = point || new THREE.Vector3(0, 0, 0);
    controls.target.copy(point);
    controls.update();
}

// ============================================================================
// EXPORT FOR USE IN OTHER SCRIPTS
// ============================================================================

// If using module system, export these:
if (typeof module !== 'undefined' && module.exports) {
    module.exports = {
        scene: scene,
        camera: camera,
        renderer: renderer,
        controls: controls,
        atoms: atoms,
        bondsArray: bondsArray,
        initializeViewer: initializeViewer,
        clearScene: clearScene,
        resetView: resetView,
        centerOnPoint: centerOnPoint,
        activeFunction: activeFunction,
        clickFunction: clickFunction,
        firstIntersectedObject: firstIntersectedObject
    };
}
