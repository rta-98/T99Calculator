function toggleGrid(homePath = {

    '/' : 'grid', 
    '/grid': '/', 
//    '/filter-grid': '/'

  }) {

      const dir1 = window.location.pathname;
      window.location.href = homePath[dir1] || '/';  

}

function toggleFilterGrid(homePath = {

    '/' : 'filter-grid', 
    '/filter-grid': '/', 
//    '/filter-grid': '/'

  }) {

      const dir2 = window.location.pathname;
      window.location.href = homePath[dir2] || '/';  

}

export  { toggleGrid, toggleFilterGrid };
