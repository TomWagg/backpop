document.addEventListener("DOMContentLoaded", function() {
    // if the secondary sidebar is empty then expand the main area and hide the sidebar
    let secondarySidebar = document.querySelector('.md-sidebar--secondary');
    if (secondarySidebar && secondarySidebar.innerText.trim() === '') {
        let mainContent = document.querySelector('.md-content');
        if (mainContent) {
            mainContent.style.marginRight = '0';
        }
        secondarySidebar.style.display = 'none';
    }

    // hide the Home title
    const h1 = document.querySelector('h1');
    if (h1 && h1.innerText.trim().startsWith('Home')) {
        h1.style.display = 'none';
    }
});