// modules glossary:
// utilities.js - shared helper functions
// smiles-api.js - handles ajax requests, submission, verification
// mol-display.js - displays single mol image
// gallery.js - displays gallery of mol images; molecular info as well
// input-form.js - smiles input box and submit button 
// display.js - 3d gui molecule visualization  

// exporting rules: 
// only export var, let, const, and functions -- only export top level. the best way to approach this is to export functions at the end of a file like so: 
// 
//
// export { name, draw, reportarea, reportperimeter }; 
//
// importing rules:  
// use the import statement followed by a comma-separated list of the features you want to import wrapped in curly braces, followed by the keyword from, followed by a module specifier. use dot syntax to signify current location. 
// import { name, draw, reportarea, reportperimeter } from "./modules/square.js"; 
//
// -- once imported ...
// you can use them just like they were defined inside the same file. for example inside of main.js: 
//
// const myCanvas = create("myCanvas", document.body, 480, 320);
// const reportList = createReportList(myCanvas.id);
// const square = draw(myCanvas.ctx, 60, 60 ,60, "blue");
// reportArea(square.length, reportList);
// reportPerim(square.length, reportList); 
//
// -- Applying to HTML ...
// specifically main.js, which is very similar to how one applies script minus the following:
// 1. Include type="module" in the <script> element, 
// <script type="module" src="main.js"></script>
// 2. Can be embedded directly into the html
// 3. BUT, import and export statements must have type="module" specified. 
 
import { showMol } from './modules/mol-display.js';
import { canonIn } from './modules/smiles-api.js';
import { toggleGrid, toggleFilterGrid } from './modules/toggle.js'; 
import { addClass, hasClass, removeClass, toggleClass, clearTheAllBtn, setTheAllBtn, clearActives, getActives, filterPNG } from './modules/gallery.js';    
import Basic3DViewer from './modules/basic-viewer.js';

function getInitialSmiles() {
  const container = document.getElementById('model3d_container') ||
    document.getElementById('model3d');
  if (container && container.dataset && container.dataset.smiles) {
    return container.dataset.smiles.trim();
  }

  if (typeof window !== 'undefined') {
    if (window.__SMILES__) {
      return String(window.__SMILES__).trim();
    }
    if (window.SMILES) {
      return String(window.SMILES).trim();
    }
  }

  const input = document.getElementById('smilesString');
  if (input && input.value) {
    return input.value.trim();
  }

  return '';
}

function initViewer() {
  const container = document.getElementById('model3d');
  if (!container) {
    return;
  }

  Basic3DViewer.initializeViewer({
    containerId: 'model3d',
    containerOuterId: 'model3d_container'
  });

  const smiles = getInitialSmiles();
  if (smiles) {
    Basic3DViewer.loadSmiles(smiles)
      .catch((err) => {
        console.error('Failed to load SMILES:', err);
      });
  }
}

document.addEventListener('DOMContentLoaded', initViewer);

export function loadSmilesString(smiles) {
  if (!smiles) {
    return Promise.resolve();
  }
  return Basic3DViewer.loadSmiles(String(smiles));
}
window.ajaxMol = function() {

  const smiIn = document.getElementById('smiles-input').value;

  showMol(smiIn);
  canonIn(smiIn);  

}

let filterContainer;
let filterBtns;

document.addEventListener("DOMContentLoaded", () => {
  filterContainer = document.getElementById("filterChildren");
  if (!filterContainer) return;

  filterBtns = filterContainer.getElementsByClassName("filterBtn");

  filterPNG("all");

  for (let i = 0; i < filterBtns.length; i++) {
    filterBtns[i].addEventListener("click", function(){
      let filter = this.getAttribute("data-filter") || "all";
      if (filter === "all") {
        clearActives();
        addClass(this, "active");
      } else {
        clearTheAllBtn();
        toggleClass(this, "active");
      }
      filterPNG();
    });
  }
});

document.addEventListener('DOMContentLoaded', () => {
  const grid = document.querySelector('.all-pngs');
  if (!grid) return;
  grid.addEventListener('click', (e) => {
    const tile = e.target.closest('[data-pdb], [data-smiles]');
    if (!tile) return;
    const pdbUrl = tile.dataset.pdb;
    const smiles = (tile.dataset.smiles || '').trim();
    if (pdbUrl) {
      Basic3DViewer.loadPDBFromUrl(pdbUrl).catch(console.error);
    } else if (smiles) {
      loadSmilesString(smiles).catch(console.error);
    }
    const viewer = document.getElementById('model3d_container') ||
document.getElementById('model3d');
    if (viewer) viewer.scrollIntoView({ behavior: 'smooth' });
  });
});

document.getElementById('toggle-btn').addEventListener('click', () => toggleGrid()); 
document.getElementById('toggle-btn-filter').addEventListener('click', () => toggleFilterGrid());
document.getElementById('toggle-btn-display')?.addEventListener('click', () => {
  const viewer = document.getElementById('model3d_container');
  if (!viewer) return;
  viewer.classList.toggle('viewer-hidden');
});


