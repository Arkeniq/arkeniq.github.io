---
title:  "🧾 Making a blog on Github Pages using Jekyll"
date: 2024-01-03 13:37:00 +0200
last_modified_at: 2024-01-25 23:26:00 +0200
description: "Getting a website up and running"
tags: [Coding, Tools]
---

Decided I should start a personal journal to keep track of my workflows and small little projects. Haven't got much now but it'll hopefully add up. I'll try to include interesting data, tools or sequencing/laboratory related information sources. I'll also use this place as my own project journal - helping me keep track of things I've learned and worked on over the years. 

Now let's talk about how this webpage came to life. The whole concept is rather easy to understand.

### Step 1: Jekyll local installation
Jekyll is a static site generator written in Ruby that works by converting plain text. You can use Markdown, HTML, CSS etc. To preview the site before publishing it, I installed **Jekyll** on my machine. You'll need Ruby as prereq. The official [Jekyll installation guide](https://jekyllrb.com/docs/installation/) has all the info.

### Step 2: Picking a theme
Look on Github/Google for Jekyll themes. I chose **Chirpy** after I found it on GitHub, and it was very straightforward to use. There's a chirpy-starter public template right [here](https://github.com/cotes2020/chirpy-starter) that you can use - just click on "Use this template" and create a new repo named `<username>.github.io.` It’s all in their documentation: [Jekyll Chirpy Theme Setup](https://chirpy.cotes.page/posts/getting-started/).

![VSCode](assets/images/VSCode.png){: w="300" .right}

### Step 3: Clone the repo
Using git, you can clone your repo locally, open it using VS Code and start making changes. Don't forget to edit `_config.yml`. Posting is easy enough, mostly using Markdown with a sprinkle of HTML and some CSS editing if you want to alter the theme. When you're done, push the changes to your Github branch and you're golden. 🥇


### Step 4: Getting a domain 
You can just use `<username>.github.io.` as your address but why not get a super cheap top level domain (TLD) for your website? I found this super cool price comparison tool called [TLD-list.com](https://tld-list.com/) where you can check domain prices and sort them by popularity or cheapest registration/renewal/transfer cost. If you want a complete list of available TLDs there's one made by the Internet Assigned Numbers Authority (IANA) linked [here](https://data.iana.org/TLD/tlds-alpha-by-domain.txt). Yes, apparently .pizza is a TLD and yes, it costs more than a pizza 🍕 to renew every year, compared to .com which is currently super cheap, at just $10/yr.

---

I'm relatively new to the Github world but having forced myself to use it to make this site has been a blast. Looking forward to using it more and incorporating more projects into it. 🙂