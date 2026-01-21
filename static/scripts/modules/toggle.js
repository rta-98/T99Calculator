function toggle(homePath = '/', gridPath = '/grid') {
    const dir = window.location.pathname;
    window.location.href = dir === homePath ? gridPath : homePath; 
}

export  { toggle };
