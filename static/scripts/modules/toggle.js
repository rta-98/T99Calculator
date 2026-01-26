function toggle(homePath = {

    '/' : 'grid', 
    '/grid': '/filter-grid', 
    '/filter-grid': '/'

  }) {

      const dir = window.location.pathname;
      window.location.href = homePath[dir] || '/';  

}

export  { toggle };
