function canonIn(smiIn) {
 const usrStr = { smi_in: smiIn };
 const canonDisRef = document.getElementById('canon-smi');
 const inputBox = document.getElementById('smiles-input');
 const errorMsg = document.getElementById('error-message')

 return fetch('/canonical-smiles', {
   method: 'POST',
   body: JSON.stringify(usrStr),
   headers: {
     'Content-Type': 'application/json' 
   }
 })
 .then(resp => resp.json())

 .then(data => {
   if (data.result) {

     canonDisRef.textContent = data.result;
     inputBox.style.backgroundColor = 'green'; 
//       inputBox.classList.remove('result'); 
//       inputBox.classList.add('error'); 
     errorMsg.textContent = '';

  } else {

     canonDisRef.textContent = '';
     inputBox.style.backgroundColor = 'red'; 
//       inputBox.classList.remove('error'); 
//       inputBox.classList.add('result'); 
     errorMsg.textContent = 'Please try again.';

  }})

   .catch(error => {
      inputBox.classList.remove('error'); 
      inputBox.classList.add('result'); 
      errorMsg.textContent  = 'Network Error';
      canonDisRef.textContent = '';
  });
}

export { canonIn }; 

