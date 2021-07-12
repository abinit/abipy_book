
# Documentation Howto 

To insert a link to an Abinit variable, use:

    {{ecut}}

or:

    {{asr_anaddb}}

if you want to refer to the asr variable of the anaddb code 
(note that in abimkdocs one uses asr@anaddb but, unfortunately this syntax is not supported my MyST).

There are several other replacements defined in _config.yml, mainly for links 
to the abipy documentation, python libraries and other useful websites.
Also in this case, the syntax is:

    {{replacement}}

The snippets directory contains (well) snippets that can included in the md file with the syntax:

    ```{include} snippets/abicheck_warning.md
    ```

Note that the path of the include file is **relative to the document**.

