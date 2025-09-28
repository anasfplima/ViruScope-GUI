// Homepage sliders
document.addEventListener("DOMContentLoaded", () => {
    const homeSection = document.querySelector(".home-text");
    if (homeSection) {
        homeSection.classList.add("visible");
    }
});
// Scroll-triggered animations
document.addEventListener('scroll', function() {
    const sections = document.querySelectorAll('#home-page .home-more');
    sections.forEach(section => {
        const rect = section.getBoundingClientRect();
        if (rect.top < window.innerHeight * 0.8) {
            section.classList.add('visible');
        }
    });
});

// Carousel control logic

document.querySelectorAll(".action-button").forEach(btn => {
  btn.addEventListener("click", function(e) {
    e.preventDefault(); // prevent default <a> behavior
    const id = e.target.id;

    switch (id) {
      case "arolitcard":
        Shiny.setInputValue("main_tab", "arolit", {priority: "event"});
        break;
      case "isopcard":
        Shiny.setInputValue("main_tab", "isop", {priority: "event"});
        break;
      case "comboscard":
        Shiny.setInputValue("main_tab", "primer_combinations", {priority: "event"});
        break;
      case "utilitiescard":
        Shiny.setInputValue("main_tab", "utilities", {priority: "event"});
        break;
    }
  });
});

document.addEventListener("DOMContentLoaded", () => {
  const radios = document.querySelectorAll(".carouselcontainer input[type='radio']");
  const labelElement = document.getElementById("labeltext");

  // Independent labels for each radio (in order)
  const customLabels = ["Articles", "Genomic data", "Primers", "Unprocessed data"];

  // Find the initially checked radio
  let current = Array.from(radios).findIndex(r => r.checked);

  function updateCarousel(index) {
    current = index;
    radios[current].checked = true;
    labelElement.textContent = customLabels[current];
    // // Re-trigger animation
    // labelElement.classList.add("animate");
    // void labelElement.offsetWidth; // forces reflow
    // labelElement.classList.remove("animate");
  }

//   function updateCarousel(index) {
//     current = index;
//     radios[current].checked = true;
//     labelElement.textContent = customLabels[current];
//     // Re-trigger animation
//     let text = document.querySelector('carousellabel animate')
//     text.classList.remove("animate");
//     void text.offsetWidth;
//     setTimeout(() => {
//       text.classList.add("animate");
//     }, 10)
//   }

  // Update when clicking the arrows
  document.querySelector(".leftarrow").addEventListener("click", () => {
    updateCarousel((current - 1 + radios.length) % radios.length); // wraps around
  });

  document.querySelector(".rightarrow").addEventListener("click", () => {
    updateCarousel((current + 1) % radios.length);
  });

  // Update when selecting a radio button directly
  radios.forEach((radio, index) => {
    radio.addEventListener("change", () => updateCarousel(index));
  });

  // Sync label on initial load
  updateCarousel(current);
});
