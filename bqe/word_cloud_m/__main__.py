from word_cloud_m import word_cloud as wc

if __name__ == '__main__':
    test_text = "Lorem Ipsum is simply dummy text of the printing and typesetting industry. Lorem Ipsum has been the industry's standard dummy text ever since the 1500s, when an unknown printer took a galley of type and scrambled it to make a type specimen book. It has survived not only five centuries, but also the leap into electronic typesetting, remaining essentially unchanged. It was popularised in the 1960s with the release of Letraset sheets containing Lorem Ipsum passages, and more recently with desktop publishing software like Aldus PageMaker including versions of Lorem Ipsum."
    outfile = 'test_out.jpg'
    strings_exclude = ['lorem','ipsum']
    wc.mk_wordcloud(words = test_text, filename_out = outfile, strings_exclude = strings_exclude)

