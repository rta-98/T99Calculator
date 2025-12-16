// Set viewing window dimensions 
var modelSize = {
    aspect: 2,
    width: 800,
    height: 400
};

var scene = new THREE.Scene();
var camera = new THREE.PerspectiveCamera(30, modelSize.aspect, .5, 200);

// Check for WebGL support
// No browser support? Use CanvasRenderer
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

var renderer;
if (webglAvailable()) {
    renderer = new THREE.WebGLRenderer({
        preserveDrawingBuffer: true,
        alpha: true
    });
} else {
    renderer = new THREE.CanvasRenderer();
}

//  Render 2D text overlays
var labelRenderer = new THREE.CSS2DRenderer();
labelRenderer.domElement.style.position = 'absolute';
labelRenderer.domElement.style.top = '0';
labelRenderer.domElement.style.pointerEvents = 'none';
labelRenderer.domElement.id = 'labelPlace';
renderer.setSize(modelSize.width, modelSize.height);
renderer.setClearColor(0xffffff, 0);

// Controls
var controls = new THREE.TrackballControls(camera, renderer.domElement);
controls.target.set(0, 0, 0);
controls.position0.set(0, 0, 4); 
controls.rotateSpeed = 3.0;
controls.zoomSpeed = 1.2;
controls.panSpeed = 0.8;
controls.noZoom = false;
controls.noPan = false;
controls.staticMoving = false;
controls.dynamicDampingFactor = 0.15;
controls.reset(); 
controls.enabled = false; 

// Lighting
var light = new THREE.DirectionalLight(0xffffff, 1);
scene.add(light);
var light_update = function () {
    'use strict';
    light.position.copy(camera.position);
};

light_update();
controls.addEventListener('change', light_update);
var mouse = new THREE.Vector2();
var raycaster = new THREE.Raycaster();
var atoms = [];
var bondsArray = [];

// Obtains first object from click event via @param 
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

var takenAction = function () {};

function activeFunction(func) {
    takenAction = function (param) {
        func(param);
    };
}

function onMouseClick(event) {
    var intersected = firstIntersectedObject(event);
    if (intersected != null) {
        takenAction(intersected.object);
    }
}

var moleculeClick = onMouseClick;

function clickFunction(func) {
    var container = document.getElementById('model3d_container');
    container.removeEventListener('mousedown', moleculeClick, true);
    moleculeClick = func;
    container.addEventListener('mousedown', moleculeClick, true);
}


var render = function () {
    'use strict';
    requestAnimationFrame(render);
    controls.update();
    renderer.render(scene, camera);
    labelRenderer.render(scene, camera);
};

render();

function clearScene() {
    for (var i = 0; i < atoms.length; i++) {
        scene.remove(atoms[i]);
    }
    
    for (var j = 0; j < bondsArray.length; j++) {
        scene.remove(bondsArray[j]);
    }
    
    while (scene.getObjectByName("arrow")) {
        var arrow = scene.getObjectByName('arrow');
        scene.remove(arrow);
    }
    
    while (scene.getObjectByName("orbital")) {
        var orb = scene.getObjectByName('orbital');
        scene.remove(orb);
    }
    
    while (scene.getObjectByName("cube")) {
        var cube = scene.getObjectByName('cube');
        scene.remove(cube);
    }
    bondsArray = [];
    atoms = [];
}

function clearSelection() {
}

// Handles canvas resizing 
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

// Initialize 3D viewer 
function initializeViewer() {
    var container = document.getElementById('model3d');
    if (!container) {
        console.error('Container element with id "model3d" not found');
        return;
    }
    
    container.appendChild(renderer.domElement);
    container.appendChild(labelRenderer.domElement);
    var canvas = renderer.domElement;
    canvas.setAttribute("id", "mol3dContext");
    canvas.style.backgroundColor = "#000000";
    var containerElement = document.getElementById('model3d_container');
    if (containerElement) {
        containerElement.addEventListener('mousedown', moleculeClick, true);
    }
    
    if (container) {
        container.addEventListener('click', function() {
            controls.enabled = true;
        });
    }
    
    window.addEventListener('resize', function() {
        var container = document.getElementById('model3d');
        if (container) {
            resize_canvas(parseInt(container.offsetWidth));
        }
    });
    
    resize_canvas(parseInt(container.offsetWidth));
    console.log('3D Viewer initialized');
}

function resetView() {
    controls.reset();
}

function centerOnPoint(point) {
    point = point || new THREE.Vector3(0, 0, 0);
    controls.target.copy(point);
    controls.update();
}

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
