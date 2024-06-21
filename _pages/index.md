---
layout: page
title: Home
id: home
permalink: /
---
<p style="padding: 3em 1em; background: #f5f7ff; border-radius: 4px;">
Hi! My name is Ryan Lee, and I am currently an undergraduate at MIT. I am planning to major in physics and considering a double major in either math or computer science.
<br> <br>
My primary interest is theoretical condensed matter physics, especially a field-theoretical approach. More specifically, I am interested in the topological phases of matter and studying them using TQFTs and CFTs.
<br> <br>
Before entering MIT, I participated 53rd International Physics Olympiad as a representative of South Korea (2024). I am also alumni of the 9th Hanseong Sonjaehan Scholarship (2021-2023).
</p>

<strong>Research</strong>
- [First research](https://rhlee430.github.io/first-research)

<strong>Recently updated notes</strong>

<ul>
  {% assign recent_notes = site.notes | sort: "last_modified_at_timestamp" | reverse %}
  {% for note in recent_notes limit: 3 %}
    <li>
      {{ note.last_modified_at | date: "%Y-%m-%d" }} — <a class="internal-link" href="{{ site.baseurl }}{{ note.url }}">{{ note.title }}</a>
    </li>
  {% endfor %}
</ul>

<style>
  .wrapper {
    max-width: 46em;
  }
</style>
