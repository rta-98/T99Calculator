/*
 * Basic 3D viewer module
 * Extracted from perry3d.js + molecule.js + filereader.js + elements.js
 *
 * Dependencies (ES module builds, loaded via CDN):
 * - https://unpkg.com/three/build/three.module.js
 * - https://unpkg.com/three/examples/jsm/controls/TrackballControls.js
 * - https://unpkg.com/three/examples/jsm/renderers/CSS2DRenderer.js
 */

import * as THREE from 'https://esm.sh/three@0.152.2';
import { TrackballControls } from
'https://esm.sh/three@0.152.2/examples/jsm/controls/TrackballControls.js';
import { CSS2DRenderer, CSS2DObject } from
'https://esm.sh/three@0.152.2/examples/jsm/renderers/CSS2DRenderer.js';

// ---------------------------------------------------------------------------
// VIEWER STATE
// ---------------------------------------------------------------------------

const modelSize = {
  aspect: 2,
  width: 800,
  height: 400
};

let scene = null;
let camera = null;
let renderer = null;
let labelRenderer = null;
let controls = null;
let light = null;
let atoms = [];
let bondsArray = [];
let mouse = null;
let raycaster = null;
let renderId = null;
let viewerContainer = null;
let viewerContainerOuter = null;
let moleculeClick = null;

// ---------------------------------------------------------------------------
// VIEWER SETUP
// ---------------------------------------------------------------------------

function webglAvailable() {
  try {
    const canvas = document.createElement('canvas');
    return !!(window.WebGLRenderingContext && (
      canvas.getContext('webgl') ||
      canvas.getContext('experimental-webgl')
    ));
  } catch (e) {
    return false;
  }
}

function initializeViewer(options = {}) {
  if (!webglAvailable()) {
    throw new Error('WebGL is not available in this environment.');
  }

  const containerId = options.containerId || 'model3d';
  const containerOuterId = options.containerOuterId || 'model3d_container';
  const container = document.getElementById(containerId);

  if (!container) {
    throw new Error(`Container element with id "${containerId}" not found`);
  }

  // Reset prior viewer if re-initializing.
  destroyViewer();

  viewerContainer = container;
  viewerContainerOuter = document.getElementById(containerOuterId) || container;

  scene = new THREE.Scene();
  camera = new THREE.PerspectiveCamera(30, modelSize.aspect, 0.5, 200);

  renderer = new THREE.WebGLRenderer({
    preserveDrawingBuffer: true,
    alpha: true
  });
  renderer.setPixelRatio(window.devicePixelRatio || 1);

  labelRenderer = new CSS2DRenderer();
  labelRenderer.domElement.style.position = 'absolute';
  labelRenderer.domElement.style.top = '0';
  labelRenderer.domElement.style.pointerEvents = 'none';
  labelRenderer.domElement.id = options.labelRendererId || 'labelPlace';

  const width = options.width || container.offsetWidth || modelSize.width;
  const height = options.height || Math.round(width / modelSize.aspect);
  modelSize.width = width;
  modelSize.height = height;
  modelSize.aspect = width / height;

  renderer.setSize(modelSize.width, modelSize.height);
  renderer.setClearColor(0xffffff, 0);
  labelRenderer.setSize(modelSize.width, modelSize.height);

  container.innerHTML = '';
  container.appendChild(renderer.domElement);
  container.appendChild(labelRenderer.domElement);

  const canvas = renderer.domElement;
  canvas.setAttribute('id', options.canvasId || 'mol3dContext');
  canvas.style.backgroundColor = '#000000';

  controls = new TrackballControls(camera, renderer.domElement);
  controls.target.set(0, 0, 0);
  controls.position0.set(0, 0, 17);
  controls.rotateSpeed = 3.0;
  controls.zoomSpeed = 1.2;
  controls.panSpeed = 0.8;
  controls.noZoom = false;
  controls.noPan = false;
  controls.staticMoving = false;
  controls.dynamicDampingFactor = 0.15;
  controls.reset();
  controls.enabled = !!options.controlsEnabled;

  light = new THREE.DirectionalLight(0xffffff, 1);
  scene.add(light);

  const lightUpdate = () => {
    light.position.copy(camera.position);
  };

  lightUpdate();
  controls.addEventListener('change', lightUpdate);

  mouse = new THREE.Vector2();
  raycaster = new THREE.Raycaster();

  moleculeClick = onMouseClick;
  if (viewerContainerOuter) {
    viewerContainerOuter.addEventListener('mousedown', moleculeClick, true);
    viewerContainerOuter.addEventListener('click', () => {
      controls.enabled = true;
    });
  }

  window.addEventListener('resize', onWindowResize);

  render();
  return true;
}

function destroyViewer() {
  if (renderId) {
    cancelAnimationFrame(renderId);
    renderId = null;
  }

  if (viewerContainerOuter && moleculeClick) {
    viewerContainerOuter.removeEventListener('mousedown', moleculeClick, true);
  }
  window.removeEventListener('resize', onWindowResize);

  if (viewerContainer) {
    viewerContainer.innerHTML = '';
  }

  scene = null;
  camera = null;
  renderer = null;
  labelRenderer = null;
  controls = null;
  light = null;
  atoms = [];
  bondsArray = [];
  mouse = null;
  raycaster = null;
  viewerContainer = null;
  viewerContainerOuter = null;
  moleculeClick = null;
}

function onWindowResize() {
  if (!viewerContainer) {
    return;
  }
  resizeCanvas(viewerContainer.offsetWidth);
}

function resizeCanvas(newWidth, newHeight) {
  if (!renderer || !camera) {
    return;
  }
  const width = newWidth || modelSize.width;
  const height = newHeight || Math.round(width / modelSize.aspect);
  modelSize.width = width;
  modelSize.height = height;
  modelSize.aspect = width / height;
  renderer.setSize(modelSize.width, modelSize.height);
  labelRenderer.setSize(modelSize.width, modelSize.height);
  camera.aspect = modelSize.aspect;
  camera.updateProjectionMatrix();
}

function render() {
  renderId = requestAnimationFrame(render);
  if (controls) {
    controls.update();
  }
  if (renderer && scene && camera) {
    renderer.render(scene, camera);
  }
  if (labelRenderer) {
    labelRenderer.render(scene, camera);
  }
}

function resetView() {
  if (controls) {
    controls.reset();
  }
}

function centerOnPoint(point) {
  if (!controls) {
    return;
  }
  const target = point || new THREE.Vector3(0, 0, 0);
  controls.target.copy(target);
  controls.update();
}

// ---------------------------------------------------------------------------
// RAYCASTING AND INTERACTION
// ---------------------------------------------------------------------------

let takenAction = () => {};

function activeFunction(func) {
  takenAction = (param) => {
    func(param);
  };
}

function clickFunction(func) {
  if (!viewerContainerOuter) {
    return;
  }
  viewerContainerOuter.removeEventListener('mousedown', moleculeClick, true);
  moleculeClick = func;
  viewerContainerOuter.addEventListener('mousedown', moleculeClick, true);
}

function onMouseClick(event) {
  const intersected = firstIntersectedObject(event);
  if (intersected != null) {
    takenAction(intersected.object);
  }
}

function firstIntersectedObject(event) {
  if (!renderer || !camera || !raycaster || !mouse) {
    return null;
  }
  event.preventDefault();

  const rect = renderer.domElement.getBoundingClientRect();
  mouse.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
  mouse.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;

  raycaster.setFromCamera(mouse, camera);

  const intersects = raycaster.intersectObjects(atoms);
  const bondIntersects = raycaster.intersectObjects(bondsArray);

  if (intersects.length > 0) {
    if (bondIntersects.length > 0) {
      const total = bondIntersects.concat(intersects);
      const tot = total.map((hit) => hit.object);
      const all = raycaster.intersectObjects(tot);
      return all[0];
    }
    return intersects[0];
  }

  if (bondIntersects.length > 0) {
    return bondIntersects[0];
  }

  return null;
}

// ---------------------------------------------------------------------------
// MOLECULE DATA STRUCTURES
// ---------------------------------------------------------------------------

function parameters() {
  if (typeof parameters.values === 'undefined') {
    parameters.values = new ParamObject();
  }
  return parameters.values;
}

class ParamObject {
  constructor() {
    this.atomSize = 0.5;
    this.atomScale = 1.0;
    this.atomResolution = 30;
    this.atomTrans = false;
    this.atomOpacity = 0.85;
    this.bondWidth = 0.08;
    this.bondResolution = 24;
    this.bondColor = 'rgb(180,180,180)';
    this.bondsOpacity = 1;
    this.bondsTrans = false;
    this.fontSize = 14;
    this.showLabels = false;
  }
}

function Mol(value) {
  const param = value || 0;
  if (typeof Mol.molecule === 'undefined') {
    Mol.molecule = 1;
  }
  if (param > 0) {
    Mol.molecule = param;
  }
  const id = Mol.molecule;

  if (typeof Mol.moly === 'undefined') {
    Mol.moly = [];
  }
  if (typeof Mol.moly[id] === 'undefined') {
    Mol.moly[id] = [];
    Mol.moly[id][0] = new MolObject();
    Mol.moly[id][0].molIndex = id;
    Mol.moly[id][0].numatoms = 0;
  }

  return Mol.moly[id];
}

class AtomObject {
  constructor() {
    this.atomicnumber = 0;
    this.x = 0.0;
    this.y = 0.0;
    this.z = 0.0;
    this.charge = 0.0;
    this.highlite = 0;
  }
}

class MolObject {
  constructor() {
    const params = parameters();
    this.molIndex = 1;
    this.numatoms = 0;
    this.AtomScale = params.atomScale;
    this.showlabels = 0;
    this.showcharges = 0;
    this.gradients = 1;
    this.highlite = 0;
    this.hide = 0;
    this.formula = '';
    this.weight = 0.0;
    this.charge = 999;
    this.center = [0, 0, 0, 0, 0];
    this.bonds = [];
  }
}

function delMolecule(value = 0) {
  const molecule = Mol(value);
  molecule[0] = new MolObject();
  molecule[0].molIndex = 1;
  molecule[0].numatoms = 0;
  for (let i = molecule.length; i > 1; i--) {
    molecule.pop();
  }
}

function addAtom(AtomicNum, x, y, z, value = 0) {
  const molecule = Mol(value);
  const bonds = molecule[0].bonds;
  molecule[0].numatoms++;
  const numatoms = molecule[0].numatoms;
  molecule[numatoms] = new AtomObject();
  molecule[numatoms].atomicnumber = AtomicNum;
  molecule[numatoms].x = x;
  molecule[numatoms].y = y;
  molecule[numatoms].z = z;
  molecule[numatoms].charge = 0.0;
  molecule[numatoms].highlite = 0;

  if (typeof bonds[numatoms] === 'undefined') {
    bonds[numatoms] = [];
  }
}

function addBond(atom1, atom2, molecule = Mol(0)) {
  const bonds = molecule[0].bonds;
  if (!bonds[atom1]) {
    bonds[atom1] = [];
  }
  if (!bonds[atom2]) {
    bonds[atom2] = [];
  }
  if (bonds[atom1].indexOf(atom2) === -1) {
    bonds[atom1].push(atom2);
  }
}

// ---------------------------------------------------------------------------
// DRAWING
// ---------------------------------------------------------------------------

function clearScene() {
  if (!scene) {
    return;
  }
  atoms.forEach((obj) => scene.remove(obj));
  bondsArray.forEach((obj) => scene.remove(obj));
  while (scene.getObjectByName('label')) {
    scene.remove(scene.getObjectByName('label'));
  }
  while (scene.getObjectByName('extra')) {
    scene.remove(scene.getObjectByName('extra'));
  }
  atoms = [];
  bondsArray = [];
}

function clearSelection() {
  // Intentionally minimal for standalone use.
}

function drawMolecule(value = 0) {
  const molecule = Mol(value);
  const params = parameters();

  if (!scene || !renderer || !camera) {
    return;
  }

  clearScene();

  const numatoms = molecule[0].numatoms;
  for (let i = 0; i < numatoms; i++) {
    const atomIndex = i + 1;
    drawAtom(atomIndex, params.atomSize, value);
  }
  bondBonded(value);
}

function drawAtom(atomNum, AtomSize, value = 0) {
  const molecule = Mol(value);
  const params = parameters();
  const AtomScale = params.atomScale;
  const A = molecule[atomNum].atomicnumber;
  const x = AtomScale * molecule[atomNum].x;
  const y = AtomScale * molecule[atomNum].y;
  const z = AtomScale * molecule[atomNum].z;
  const r = AtomScale * AtomSize * element(A, 'radius');

  const geometry = new THREE.SphereGeometry(r, params.atomResolution, params.atomResolution);
  const material = new THREE.MeshPhongMaterial({
    color: element(A, 'color'),
    transparent: params.atomTrans,
    opacity: params.atomOpacity
  });
  const atom = new THREE.Mesh(geometry, material);
  atom.position.set(x, y, z);
  atom.name = atomNum;
  scene.add(atom);
  atoms.push(atom);

  if (params.showLabels) {
    atomLabel(A, atomNum, value);
  }
}

function drawBond(atom1, atom2, value = 0) {
  const params = parameters();
  const molecule = Mol(value);
  const AtomScale = params.atomScale;
  const BondWidth = params.bondWidth;
  const resolution = params.bondResolution;
  const BondColor = params.bondColor;
  const opac = params.bondsOpacity;
  const trans = params.bondsTrans;

  const x1 = AtomScale * molecule[atom1].x;
  const y1 = AtomScale * molecule[atom1].y;
  const z1 = AtomScale * molecule[atom1].z;
  const x2 = AtomScale * molecule[atom2].x;
  const y2 = AtomScale * molecule[atom2].y;
  const z2 = AtomScale * molecule[atom2].z;

  const point1 = new THREE.Vector3(x1, y1, z1);
  const point2 = new THREE.Vector3(x2, y2, z2);
  const dist = point1.distanceTo(point2);
  const bondGeo = new THREE.CylinderGeometry(BondWidth, BondWidth, dist, resolution);
  const bondMat = new THREE.MeshPhongMaterial({
    color: BondColor,
    transparent: trans,
    opacity: opac
  });
  const bondMesh = cylinderMesh(point1, point2, bondMat, bondGeo);
  bondMesh.name = `bond_${atom1}-${atom2}`;
  bondsArray.push(bondMesh);
  scene.add(bondMesh);
}

function bondBonded(value = 0) {
  const molecule = Mol(value);
  const bonds = molecule[0].bonds;
  const numatoms = molecule[0].numatoms;

  for (let i = 0; i < numatoms; i++) {
    const atomIndex = i + 1;
    if (!bonds[atomIndex]) {
      continue;
    }
    for (let j = 0; j < bonds[atomIndex].length; j++) {
      const connection = bonds[atomIndex][j];
      if (connection > atomIndex) {
        drawBond(atomIndex, connection, value);
      }
    }
  }
}

function atomLabel(A, i, value) {
  const params = parameters();
  const molecule = Mol(value);
  const AtomScale = params.atomScale;
  const label = element(A, 'symbol');
  const color = element(A, 'color');
  const atom = molecule[i];
  const x = atom.x * AtomScale;
  const y = atom.y * AtomScale;
  const z = atom.z * AtomScale;

  const text = document.createElement('div');
  text.className = 'label';
  text.style.color = color;
  text.textContent = `${label} ${i}`;

  const cssLabel = new CSS2DObject(text);
  cssLabel.position.set(x, y, z);
  cssLabel.name = 'label';
  scene.add(cssLabel);
}

function cylinderMesh(pointX, pointY, material, edgeGeometry) {
  const orientation = new THREE.Matrix4();
  orientation.lookAt(pointX, pointY, new THREE.Object3D().up);
  orientation.multiply(new THREE.Matrix4().set(
    1, 0, 0, 0,
    0, 0, 1, 0,
    0, -1, 0, 0,
    0, 0, 0, 1
  ));
  const edge = new THREE.Mesh(edgeGeometry, material);
  edge.applyMatrix4(orientation);
  edge.position.x = (pointY.x + pointX.x) / 2;
  edge.position.y = (pointY.y + pointX.y) / 2;
  edge.position.z = (pointY.z + pointX.z) / 2;
  return edge;
}

// ---------------------------------------------------------------------------
// SMILES / PDB LOADING
// ---------------------------------------------------------------------------

function readPDBText(text, value) {
  const lines = text.split('\n');
  return readPDBLines(lines, value);
}

function readPDBLines(lines, value = 0) {
  const molecule = Mol(value);
  let connected = false;

  delMolecule(value);

  for (let j = 0; j < lines.length; j++) {
    const line = lines[j];
    if (!line) {
      continue;
    }
    const prefix = line.substr(0, 6).trim();
    if (prefix === 'ATOM' || prefix === 'HETATM') {
      let symbol = line.substr(76, 2).trim();
      if (!symbol) {
        const parts = line.trim().split(/\s+/);
        if (prefix === 'HETATM' && parts.length >= 11) {
          symbol = parts[10];
        } else if (parts.length >= 3) {
          symbol = parts[2];
        }
      }
      if (symbol && symbol.length > 1) {
        symbol = symbol.charAt(0).toUpperCase() + symbol.charAt(1).toLowerCase();
      }
      const A = lookupNumber(symbol);
      if (!A) {
        continue;
      }
      const x = parseFloat(line.substr(30, 8));
      const y = parseFloat(line.substr(38, 8));
      const z = parseFloat(line.substr(46, 8));
      addAtom(A, x, y, z, value);
    } else if (prefix === 'CONECT') {
      connected = true;
      const tokens = line.trim().split(/\s+/);
      const from = parseInt(tokens[1], 10);
      if (Number.isNaN(from)) {
        continue;
      }
      for (let k = 2; k < tokens.length; k++) {
        const to = parseInt(tokens[k], 10);
        if (!Number.isNaN(to) && to !== from) {
          addBond(from, to, molecule);
          addBond(to, from, molecule);
        }
      }
    } else if (prefix === 'END') {
      break;
    }
  }

  if (!connected) {
    const numatoms = molecule[0].numatoms;
    for (let i = 1; i < numatoms; i++) {
      for (let n = i + 1; n <= numatoms; n++) {
        const bond = 1.2 * (element(molecule[i].atomicnumber, 'radius') +
          element(molecule[n].atomicnumber, 'radius'));
        const dx = molecule[i].x - molecule[n].x;
        const dy = molecule[i].y - molecule[n].y;
        const dz = molecule[i].z - molecule[n].z;
        const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (r <= bond) {
          addBond(i, n, molecule);
          addBond(n, i, molecule);
        }
      }
    }
  }

  molecule[0].center = moleculeCenter(value);
  return molecule;
}

function loadPDBFromUrl(url, value) {
  return fetchText(url).then((text) => {
    readPDBText(text, value);
    translateMolecule('center', { redraw: true });
    return Mol(value || 0);
  });
}

function loadSmiles(smiles, value) {
  const encoded = encodeURIComponent(smiles);
  const url = `https://cactus.nci.nih.gov/ncidb2.2/nci2.2.tcl?op1=fs&method1=ens&data1=${encoded}&output=pdb&maxhits=1&nomsg=1`;

  return fetchText(url).then((text) => {
    const firstLine = text.split('\n')[0] || '';
    if (firstLine.toLowerCase().includes('html')) {
      throw new Error('SMILES lookup returned no PDB result.');
    }
    readPDBText(text, value);
    translateMolecule('center', { redraw: true });
    return Mol(value || 0);
  });
}

function fetchText(url) {
  if (typeof fetch === 'function') {
    return fetch(url).then((res) => {
      if (!res.ok) {
        throw new Error(`Request failed with status ${res.status}`);
      }
      return res.text();
    });
  }

  return new Promise((resolve, reject) => {
    const xhr = new XMLHttpRequest();
    xhr.open('GET', url, true);
    xhr.onreadystatechange = () => {
      if (xhr.readyState === XMLHttpRequest.DONE) {
        if (xhr.status >= 200 && xhr.status < 300) {
          resolve(xhr.responseText);
        } else {
          reject(new Error(`Request failed with status ${xhr.status}`));
        }
      }
    };
    xhr.send();
  });
}

// ---------------------------------------------------------------------------
// UTILITIES
// ---------------------------------------------------------------------------

function lookupNumber(symbol) {
  for (let i = 1; i < element(1, 'max'); i++) {
    if (element(i, 'symbol') === symbol) {
      return i;
    }
  }
  return 0;
}

function lookupSymbol(number) {
  if (number < element(1, 'max')) {
    return element(number, 'symbol');
  }
  return '';
}

function distanceBetween3D(x1, y1, z1, x2, y2, z2) {
  return Math.sqrt((x1 - x2) * (x1 - x2) +
    (y1 - y2) * (y1 - y2) +
    (z1 - z2) * (z1 - z2));
}

function moleculeCenter(value = 0) {
  const molecule = Mol(value);
  let x;
  let y;
  let z;
  let maxX = 0;
  let maxY = 0;
  let maxZ = 0;
  let minX = 0;
  let minY = 0;
  let minZ = 0;
  let sumX = 0;
  let sumY = 0;
  let sumZ = 0;

  for (let i = 1; i < molecule.length; i++) {
    x = parseFloat(molecule[i].x);
    y = parseFloat(molecule[i].y);
    z = parseFloat(molecule[i].z);
    if (x > maxX) { maxX = x; }
    if (x < minX) { minX = x; }
    if (y > maxY) { maxY = y; }
    if (y < minY) { minY = y; }
    if (z > maxZ) { maxZ = z; }
    if (z < minZ) { minZ = z; }
    sumX += x;
    sumY += y;
    sumZ += z;
  }

  const avgX = sumX / molecule[0].numatoms;
  const avgY = sumY / molecule[0].numatoms;
  const avgZ = sumZ / molecule[0].numatoms;
  const bigX = Math.max(Math.abs(minX), Math.abs(maxX));
  const bigY = Math.max(Math.abs(minY), Math.abs(maxY));
  const bigZ = Math.max(Math.abs(minZ), Math.abs(maxZ));

  if ((bigX === maxX) || (bigX === -maxX)) { x = maxX; }
  if ((bigX === minX) || (bigX === -minX)) { x = minX; }
  if ((bigY === maxY) || (bigY === -maxY)) { y = maxY; }
  if ((bigY === minY) || (bigY === -minY)) { y = minY; }
  if ((bigZ === maxZ) || (bigZ === -maxZ)) { z = maxZ; }
  if ((bigZ === minZ) || (bigZ === -minZ)) { z = minZ; }

  return [avgX, avgY, avgZ, x, y, z];
}

function translateMolecule(distance, options = {}) {
  const redraw = options.redraw !== false;
  const molecule = Mol();

  if (distance === 'center') {
    const center = molecule[0].center;
    for (let i = 1; i < molecule.length; i++) {
      molecule[i].x -= center[0];
      molecule[i].y -= center[1];
      molecule[i].z -= center[2];
    }
    molecule[0].center = moleculeCenter();
    if (redraw) {
      drawMolecule(0);
    }
    return;
  }

  if (typeof distance !== typeof []) {
    return;
  }
  if (typeof distance[0] === 'undefined') {
    return;
  }
  for (let j = 1; j < molecule.length; j++) {
    molecule[j].x += distance[0];
    molecule[j].y += distance[1];
    molecule[j].z += distance[2];
  }
  molecule[0].center = moleculeCenter();
  if (redraw) {
    drawMolecule(0);
  }
}

// ---------------------------------------------------------------------------
// ELEMENT DATA (from elements.js)
// ---------------------------------------------------------------------------

function element(Z, param) {
  if (typeof element.elem === 'undefined') {
    element.elem = [];
    element.elem[0] = addElement(0, 'X', 's', 0, 0, 0, 0, 0, 0, 0);
    element.elem[1] = addElement(1, 'H', 's', 1, 1.0079, 31, 255, 255, 255, 2.2);
    element.elem[2] = addElement(2, 'He', 's', 2, 4.0026, 28, 217, 255, 255, 0);
    element.elem[3] = addElement(3, 'Li', 's', 1, 6.941, 128, 204, 128, 255, 0.98);
    element.elem[4] = addElement(4, 'Be', 's', 2, 9.0122, 96, 194, 255, 0, 1.57);
    element.elem[5] = addElement(5, 'B', 'p', 3, 10.811, 84, 255, 181, 181, 2.04);
    element.elem[6] = addElement(6, 'C', 'p', 4, 12.0107, 76, 144, 144, 144, 2.55);
    element.elem[7] = addElement(7, 'N', 'p', 5, 14.0067, 71, 48, 80, 248, 3.04);
    element.elem[8] = addElement(8, 'O', 'p', 6, 15.9994, 66, 255, 13, 13, 3.44);
    element.elem[9] = addElement(9, 'F', 'p', 7, 18.9984, 57, 144, 224, 80, 3.98);
    element.elem[10] = addElement(10, 'Ne', 'p', 8, 20.1797, 58, 179, 227, 245, 0);
    element.elem[11] = addElement(11, 'Na', 's', 1, 22.9897, 166, 171, 92, 242, 0.93);
    element.elem[12] = addElement(12, 'Mg', 's', 2, 24.305, 141, 138, 255, 0, 1.31);
    element.elem[13] = addElement(13, 'Al', 'p', 3, 26.9815, 121, 191, 166, 166, 1.61);
    element.elem[14] = addElement(14, 'Si', 'p', 4, 28.0855, 111, 240, 200, 160, 1.9);
    element.elem[15] = addElement(15, 'P', 'p', 5, 30.9738, 107, 255, 128, 0, 2.19);
    element.elem[16] = addElement(16, 'S', 'p', 6, 32.065, 105, 255, 255, 48, 2.58);
    element.elem[17] = addElement(17, 'Cl', 'p', 7, 35.453, 102, 31, 240, 31, 3.16);
    element.elem[18] = addElement(18, 'Ar', 'p', 8, 39.948, 106, 128, 209, 227, 0);
    element.elem[19] = addElement(19, 'K', 's', 1, 39.0983, 203, 143, 64, 212, 0.82);
    element.elem[20] = addElement(20, 'Ca', 's', 2, 40.078, 176, 61, 255, 0, 1);
    element.elem[21] = addElement(21, 'Sc', 'd', 3, 44.9559, 170, 230, 230, 230, 1.36);
    element.elem[22] = addElement(22, 'Ti', 'd', 4, 47.867, 160, 191, 194, 199, 1.54);
    element.elem[23] = addElement(23, 'V', 'd', 5, 50.9415, 153, 166, 166, 171, 1.63);
    element.elem[24] = addElement(24, 'Cr', 'd', 6, 51.9961, 139, 138, 153, 199, 1.66);
    element.elem[25] = addElement(25, 'Mn', 'd', 7, 54.938, 139, 156, 122, 199, 1.55);
    element.elem[26] = addElement(26, 'Fe', 'd', 8, 55.845, 132, 224, 102, 51, 1.83);
    element.elem[27] = addElement(27, 'Co', 'd', 9, 58.9332, 126, 240, 144, 160, 1.88);
    element.elem[28] = addElement(28, 'Ni', 'd', 10, 58.6934, 124, 80, 208, 80, 1.91);
    element.elem[29] = addElement(29, 'Cu', 'd', 11, 63.546, 132, 200, 128, 51, 1.9);
    element.elem[30] = addElement(30, 'Zn', 'd', 2, 65.39, 122, 125, 128, 176, 1.65);
    element.elem[31] = addElement(31, 'Ga', 'p', 3, 69.723, 122, 194, 143, 143, 1.81);
    element.elem[32] = addElement(32, 'Ge', 'p', 4, 72.64, 120, 102, 143, 143, 2.01);
    element.elem[33] = addElement(33, 'As', 'p', 5, 74.9216, 119, 189, 128, 227, 2.18);
    element.elem[34] = addElement(34, 'Se', 'p', 6, 78.96, 120, 255, 161, 0, 2.55);
    element.elem[35] = addElement(35, 'Br', 'p', 7, 79.904, 120, 166, 41, 41, 2.96);
    element.elem[36] = addElement(36, 'Kr', 'p', 8, 83.8, 116, 92, 184, 209, 3);
    element.elem[37] = addElement(37, 'Rb', 's', 1, 85.4678, 220, 112, 46, 176, 0.82);
    element.elem[38] = addElement(38, 'Sr', 's', 2, 87.62, 195, 0, 255, 0, 0.95);
    element.elem[39] = addElement(39, 'Y', 'd', 3, 88.9059, 190, 148, 255, 255, 1.22);
    element.elem[40] = addElement(40, 'Zr', 'd', 4, 91.224, 175, 148, 224, 224, 1.33);
    element.elem[41] = addElement(41, 'Nb', 'd', 5, 92.9064, 164, 115, 194, 201, 1.6);
    element.elem[42] = addElement(42, 'Mo', 'd', 6, 95.94, 154, 84, 181, 181, 2.16);
    element.elem[43] = addElement(43, 'Tc', 'd', 7, 98, 147, 59, 158, 158, 1.9);
    element.elem[44] = addElement(44, 'Ru', 'd', 8, 101.07, 146, 36, 143, 143, 2.2);
    element.elem[45] = addElement(45, 'Rh', 'd', 9, 102.9055, 142, 10, 125, 140, 2.28);
    element.elem[46] = addElement(46, 'Pd', 'd', 10, 106.42, 139, 0, 105, 133, 2.2);
    element.elem[47] = addElement(47, 'Ag', 'd', 11, 107.8682, 145, 192, 192, 192, 1.93);
    element.elem[48] = addElement(48, 'Cd', 'd', 2, 112.411, 144, 255, 217, 143, 1.69);
    element.elem[49] = addElement(49, 'In', 'p', 3, 114.818, 142, 166, 117, 115, 1.78);
    element.elem[50] = addElement(50, 'Sn', 'p', 4, 118.71, 139, 102, 128, 128, 1.96);
    element.elem[51] = addElement(51, 'Sb', 'p', 5, 121.76, 139, 158, 99, 181, 2.05);
    element.elem[52] = addElement(52, 'Te', 'p', 6, 127.6, 138, 212, 122, 0, 2.1);
    element.elem[53] = addElement(53, 'I', 'p', 7, 126.9045, 139, 148, 0, 148, 2.66);
    element.elem[54] = addElement(54, 'Xe', 'p', 8, 131.293, 140, 66, 158, 176, 2.6);
    element.elem[55] = addElement(55, 'Cs', 's', 1, 132.9055, 244, 87, 23, 143, 0.79);
    element.elem[56] = addElement(56, 'Ba', 's', 2, 137.327, 215, 0, 201, 0, 0.89);
    element.elem[57] = addElement(57, 'La', 'd', 3, 138.9055, 207, 112, 212, 255, 1.1);
    element.elem[58] = addElement(58, 'Ce', 'f', 4, 140.116, 204, 255, 255, 199, 1.12);
    element.elem[59] = addElement(59, 'Pr', 'f', 5, 140.9077, 203, 217, 255, 199, 1.13);
    element.elem[60] = addElement(60, 'Nd', 'f', 6, 144.24, 201, 199, 255, 199, 1.14);
    element.elem[61] = addElement(61, 'Pm', 'f', 7, 145, 199, 163, 255, 199, 0);
    element.elem[62] = addElement(62, 'Sm', 'f', 8, 150.36, 198, 143, 255, 199, 1.17);
    element.elem[63] = addElement(63, 'Eu', 'f', 9, 151.964, 198, 97, 255, 199, 0);
    element.elem[64] = addElement(64, 'Gd', 'f', 10, 157.25, 196, 69, 255, 199, 1.2);
    element.elem[65] = addElement(65, 'Tb', 'f', 11, 158.9253, 194, 48, 255, 199, 0);
    element.elem[66] = addElement(66, 'Dy', 'f', 12, 162.5, 192, 31, 255, 199, 1.22);
    element.elem[67] = addElement(67, 'Ho', 'f', 13, 164.9303, 192, 0, 255, 156, 1.23);
    element.elem[68] = addElement(68, 'Er', 'f', 14, 167.259, 189, 0, 230, 117, 1.24);
    element.elem[69] = addElement(69, 'Tm', 'f', 15, 168.9342, 190, 0, 212, 82, 1.25);
    element.elem[70] = addElement(70, 'Yb', 'f', 16, 173.04, 187, 0, 191, 56, 0);
    element.elem[71] = addElement(71, 'Lu', 'f', 17, 174.967, 187, 0, 171, 36, 1.27);
    element.elem[72] = addElement(72, 'Hf', 'd', 4, 178.49, 175, 77, 194, 255, 1.3);
    element.elem[73] = addElement(73, 'Ta', 'd', 5, 180.9479, 170, 77, 166, 255, 1.5);
    element.elem[74] = addElement(74, 'W', 'd', 6, 183.84, 162, 33, 148, 214, 2.36);
    element.elem[75] = addElement(75, 'Re', 'd', 7, 186.207, 151, 38, 125, 171, 1.9);
    element.elem[76] = addElement(76, 'Os', 'd', 8, 190.23, 144, 38, 102, 150, 2.2);
    element.elem[77] = addElement(77, 'Ir', 'd', 9, 192.217, 141, 23, 84, 135, 2.2);
    element.elem[78] = addElement(78, 'Pt', 'd', 10, 195.078, 136, 208, 208, 224, 2.28);
    element.elem[79] = addElement(79, 'Au', 'd', 11, 196.9665, 136, 255, 209, 35, 2.54);
    element.elem[80] = addElement(80, 'Hg', 'd', 2, 200.59, 132, 184, 184, 208, 2);
    element.elem[81] = addElement(81, 'Tl', 'p', 3, 204.3833, 145, 166, 84, 77, 1.62);
    element.elem[82] = addElement(82, 'Pb', 'p', 4, 207.2, 146, 87, 89, 97, 2.33);
    element.elem[83] = addElement(83, 'Bi', 'p', 5, 208.9804, 148, 158, 79, 181, 2.02);
    element.elem[84] = addElement(84, 'Po', 'p', 6, 209, 140, 171, 92, 0, 2);
    element.elem[85] = addElement(85, 'At', 'p', 7, 210, 150, 117, 79, 69, 2.2);
    element.elem[86] = addElement(86, 'Rn', 'p', 8, 222, 150, 66, 130, 150, 0);
    element.elem[87] = addElement(87, 'Fr', 's', 1, 223, 260, 66, 0, 102, 0.7);
    element.elem[88] = addElement(88, 'Ra', 's', 2, 226, 221, 0, 125, 0, 0.9);
    element.elem[89] = addElement(89, 'Ac', 'd', 3, 227, 215, 112, 171, 250, 1.1);
    element.elem[90] = addElement(90, 'Th', 'f', 4, 232.0381, 206, 0, 186, 255, 1.3);
    element.elem[91] = addElement(91, 'Pa', 'f', 5, 231.0359, 200, 0, 161, 255, 1.5);
    element.elem[92] = addElement(92, 'U', 'f', 6, 238.0289, 196, 0, 143, 255, 1.38);
    element.elem[93] = addElement(93, 'Np', 'f', 7, 237, 190, 0, 128, 255, 1.36);
    element.elem[94] = addElement(94, 'Pu', 'f', 8, 244, 187, 0, 107, 255, 1.28);
    element.elem[95] = addElement(95, 'Am', 'f', 9, 243, 180, 84, 92, 242, 1.3);
    element.elem[96] = addElement(96, 'Cm', 'f', 10, 247, 169, 120, 92, 227, 1.3);
    element.elem[97] = addElement(97, 'Bk', 'f', 11, 247, 168, 138, 79, 227, 1.3);
    element.elem[98] = addElement(98, 'Cf', 'f', 12, 251, 168, 161, 54, 212, 1.3);
    element.elem[99] = addElement(99, 'Es', 'f', 13, 252, 165, 179, 31, 212, 1.3);
    element.elem[100] = addElement(100, 'Fm', 'f', 14, 257, 167, 179, 31, 186, 1.3);
    element.elem[101] = addElement(101, 'Md', 'f', 15, 258, 173, 179, 13, 166, 1.3);
    element.elem[102] = addElement(102, 'No', 'f', 16, 259, 176, 189, 13, 135, 1.3);
    element.elem[103] = addElement(103, 'Lr', 'f', 17, 262, 161, 199, 0, 102, 0);
    element.elem[104] = addElement(104, 'Rf', 'd', 4, 267, 157, 204, 0, 89, 0);
    element.elem[105] = addElement(105, 'Db', 'd', 5, 268, 149, 209, 0, 79, 0);
    element.elem[106] = addElement(106, 'Sg', 'd', 6, 269, 143, 217, 0, 69, 0);
    element.elem[107] = addElement(107, 'Bh', 'd', 7, 270, 141, 224, 0, 56, 0);
    element.elem[108] = addElement(108, 'Hs', 'd', 8, 269, 134, 230, 0, 46, 0);
    element.elem[109] = addElement(109, 'Mt', 'd', 9, 278, 129, 235, 0, 38, 0);
    element.elem[110] = addElement(110, 'Ds', 'd', 10, 281, 0, 0, 0, 28, 0);
    element.elem[111] = addElement(111, 'Rg', 'd', 11, 281, 0, 0, 0, 28, 0);
    element.elem[112] = addElement(112, 'Cn', 'd', 12, 285, 0, 0, 0, 28, 0);
    element.elem[113] = addElement(113, 'Uut', 'p', 3, 286, 0, 0, 0, 28, 0);
    element.elem[114] = addElement(114, 'Fl', 'p', 4, 289, 0, 0, 0, 28, 0);
    element.elem[115] = addElement(115, 'Uup', 'p', 5, 288, 0, 0, 0, 28, 0);
    element.elem[116] = addElement(116, 'Lv', 'p', 6, 293, 0, 0, 0, 28, 0);
    element.elem[117] = addElement(117, 'Uus', 'p', 7, 294, 0, 0, 0, 28, 0);
    element.elem[118] = addElement(118, 'Uuo', 'p', 8, 294, 0, 0, 0, 28, 0);
  }

  if (Z < 0 || Z > element.elem.length) {
    return 0;
  }

  switch (param) {
    case 'symbol':
      return element.elem[Z].symbol;
    case 'block':
      return element.elem[Z].block;
    case 'valence':
      return element.elem[Z].valence;
    case 'mass':
      return element.elem[Z].mass;
    case 'radius':
      return element.elem[Z].radius;
    case 'EN':
      return element.elem[Z].EN;
    case 'color':
      return element.elem[Z].color;
    case 'gradient':
      return element.elem[Z].gradient;
    case 'label':
      return element.elem[Z].label;
    case 'max':
      return element.elem.length;
    default:
      return 0;
  }
}

class ElementObject {
  constructor() {
    this.symbol = '?';
    this.block = '?';
    this.valence = 0;
    this.mass = 0.0;
    this.radius = 0.0;
    this.EN = 0.0;
    this.color = 'rgb(88,88,88)';
    this.gradient = 'rgb(88,88,88)';
    this.label = 'rgb(255,255,255)';
  }
}

function addElement(Z, Symbol, block, valence, mass, radius, red, green, blue, EN) {
  const el = new ElementObject();
  el.symbol = Symbol;
  el.block = block;
  el.valence = valence;
  el.mass = mass;
  el.radius = (radius <= 0) ? 0.1 : radius / 100.0;
  el.color = `rgb(${red},${green},${blue})`;
  el.gradient = `rgb(${Math.floor(red / 2)},${Math.floor(green / 2)},${Math.floor(blue / 2)})`;
  const sum = red + green + blue;
  el.label = (sum > 384) ? 'rgb(0,0,0)' : 'rgb(255,255,255)';
  el.EN = EN;
  return el;
}

// ---------------------------------------------------------------------------
// EXPORTS
// ---------------------------------------------------------------------------

const Basic3DViewer = {
  initializeViewer,
  destroyViewer,
  resizeCanvas,
  resetView,
  centerOnPoint,
  clearScene,
  clearSelection,
  drawMolecule,
  readPDBText,
  readPDBLines,
  loadPDBFromUrl,
  loadSmiles,
  Mol,
  parameters,
  element,
  addAtom,
  addBond,
  lookupNumber,
  lookupSymbol,
  activeFunction,
  clickFunction,
  firstIntersectedObject,
  get scene() { return scene; },
  get camera() { return camera; },
  get renderer() { return renderer; },
  get labelRenderer() { return labelRenderer; },
  get controls() { return controls; },
  get atoms() { return atoms; },
  get bondsArray() { return bondsArray; },
  get modelSize() { return modelSize; }
};

export {
  Basic3DViewer,
  initializeViewer,
  destroyViewer,
  resizeCanvas,
  resetView,
  centerOnPoint,
  clearScene,
  clearSelection,
  drawMolecule,
  readPDBText,
  readPDBLines,
  loadPDBFromUrl,
  loadSmiles,
  Mol,
  parameters,
  element,
  addAtom,
  addBond,
  lookupNumber,
  lookupSymbol,
  activeFunction,
  clickFunction,
  firstIntersectedObject,
  modelSize
};

export default Basic3DViewer;
