export function showMol(smiIn) {
  const usrStr = { smi_in: smiIn };

  fetch('/smiles-structure', { // SENDS request to server
    method: 'POST', 
    body: JSON.stringify(usrStr),
    headers: {
      'Content-Type': 'application/json' 
    }
  }) // fetch instead of XMLHttpRequest 

  .then(function(resp) { // RECIEVES reponse from server
    return resp.blob(); 
    }) 
  .then(function(blob) { // Processess response (the blob) 
    const url = URL.createObjectURL(blob);
    document.getElementById('img').src = url;
    }) // sets image source to the blob URL 
  .catch(error => console.error('Error:', error));
}


