// Modules Glossary:
// utilities.js - shared helper functions
// smiles-api.js - handles AJAX requests, submission, verification
// mol-display.js - displays single mol image
// gallery.js - displays gallery of mol images; molecular info as well
// input-form.js - SMILES input box and submit button 
// display.js - 3D Gui molecule visualization  

// Exporting Rules: 
// Only export var, let, const, and functions -- only export top level. The best way to approach this is to export functions at the end of a file like so: 
// 
//
// export { name, draw, reportArea, reportPerimeter }; 
//
// Importing Rules:  
// Use the import statement followed by a comma-separated list of the features you want to import wrapped in curly braces, followed by the keyword from, followed by a module specifier. Use dot syntax to signify current location. 
// import { name, draw, reportArea, reportPerimeter } from "./modules/square.js"; 
//
// -- Once Imported ...
// you can use them just like they were defined inside the same file. For example inside of main.js: 
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
import { toggle } from './modules/toggle.js'; 

window.ajaxMol = function() {

  const smiIn = document.getElementById('smiles-input').value;

  showMol(smiIn);
  canonIn(smiIn);  

}

document.getElementById('toggle-btn').addEventListener('click', () => toggle()); 



