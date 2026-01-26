//filterCat("all"); 
filterPNG("all"); 
//filterMode(); 

function filterPNG(c) {
  let i, png, selection; 
  if (typeof c === "string" && c == "all") {
    clearActives(); 
    setTheAllBtn(); 
  } 
  
  let pngs = document.querySelectorAll(".all-pngs .filterBtn");
  let actives = getActives() 

  for (i = 0; i < pngs.length; i++) { 
    let png = pngs[i] 
    let matches = (actives.length === 0);  

    if (actives.length > 0) {
      matches = true;
      for (let k = 0; k < actives.length; k++) {
        if (!hasClass(png, actives[k])) { 
          matches = false; 
          break; 
        }
      }
    }
    if (matches) addClass(png, "show"); // where png is the mol url's descriptor
    else removeClass(png, "show"); 
  }
}


function getActives() {
  let elements = document.getElementsByClassName("filterBtn");
  let activesArray = [];
  for (let i = 0; i < elements.length; i++) {
    if (hasClass(elements[i], "active")) {
      let filter = elements[i].getAttribute("data-filter");
      if (filter && filter !== "all") activesArray.push(filter); // adds one or more elements to array and mutates in places
    }
  }
  return activesArray;
}

function clearActives() {
  let elements = document.getElementsByClassName("filterBtn");
  for (let i = 0; i < elements.length; i++) {
    removeClass(elements[i], "active");
  }
} 

function setTheAllBtn() {
  let elements = document.getElementsByClassName("filterBtn");
  for (let i = 0; i < elements.length; i++) {
    if (elements[i].getAttribute("data-filter") === "all") {
      addClass(elements[i], "active");
      break;
    }
  }
}

function clearTheAllBtn() {
  let elements = document.getElementsByClassName("filterBtn");
  for (let i = 0; i < elements.length; i++) {
    if (elements[i].getAttribute("data-filter") === "all") {
      removeClass(elements[i], "active");
      break;
    }
  }
} 

function toggleClass(element, name) {
  if (hasClass(element, name)) removeClass(element, name);
  else addClass(element, name);
} 

function addClass(element, name) {
  let i, arr1, arr2;
  arr1 = element.className.split(" ");
  arr2 = name.split(" ");
  for (i = 0; i < arr2.length; i++) {
    if (arr1.indexOf(arr2[i]) == -1) {element.className += " " + arr2[i];}
  }

}

function removeClass(element, name) {
  let i, arr1, arr2;
  arr1 = element.className.split(" ");
  arr2 = name.split(" ");
  for (i = 0; i < arr2.length; i++) {
    while (arr1.indexOf(arr2[i]) > -1) {
      arr1.splice(arr1.indexOf(arr2[i]), 1);
    }
  }
  element.className = arr1.join(" ");
}

function hasClass(element, name) { 
  return (" " + element.className + " ").indexOf(" " + name + " ") > -1;
} 

export { addClass, hasClass, removeClass, toggleClass, clearTheAllBtn, setTheAllBtn, clearActives, getActives, filterPNG };   
