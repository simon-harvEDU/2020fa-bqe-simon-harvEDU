import pandas as pd

import datetime

from wordcloud import WordCloud, STOPWORDS 
import matplotlib.pyplot as plt 
import nltk 
import re

from sklearn.feature_extraction.text import CountVectorizer


# ideas from https://stackoverflow.com/questions/49537474/wordcloud-of-bigram-using-python

def mk_wordcloud(words, filename_out, strings_exclude):
    """
    words : string
    filename_out : string
    strings_exlucde list : ['xxx', 'yyy'] # If you want to remove any particular word form text which does not contribute much in meaning
    """
    WNL = nltk.WordNetLemmatizer()
    text = words
    # Lowercase and tokenize
    text = text.lower()
    # Remove single quote early since it causes problems with the tokenizer.
    text = text.replace("'", "")
    # Remove numbers from text
    #remove_digits = str.maketrans('', '', digits)
    #text = text.translate(remove_digits)
    tokens = nltk.word_tokenize(text)
    text1 = nltk.Text(tokens)

    # Remove extra chars and remove stop words.
    text_content = [''.join(re.split("[ .,;:!?‘’``''@#$%^_&*()<>{}~\n\t\\\-]", word)) for word in text1]

    #set the stopwords list
    stopwords_wc = set(STOPWORDS)
    
    customised_words = strings_exclude 

    new_stopwords = stopwords_wc.union(customised_words)
    text_content = [word for word in text_content if word not in new_stopwords]

    # After the punctuation above is removed it still leaves empty entries in the list.
    text_content = [s for s in text_content if len(s) != 0]

    # Best to get the lemmas of each word to reduce the number of similar words
    text_content = [WNL.lemmatize(t) for t in text_content]

    nltk_tokens = nltk.word_tokenize(text)  
    bigrams_list = list(nltk.bigrams(text_content))
    #print(bigrams_list)
    dictionary2 = [' '.join(tup) for tup in bigrams_list]
    #print (dictionary2)

    #Using count vectoriser to view the frequency of bigrams
    vectorizer = CountVectorizer(ngram_range=(2, 2))
    bag_of_words = vectorizer.fit_transform(dictionary2)
    vectorizer.vocabulary_
    sum_words = bag_of_words.sum(axis=0) 
    words_freq = [(word, sum_words[0, idx]) for word, idx in vectorizer.vocabulary_.items()]
    words_freq =sorted(words_freq, key = lambda x: x[1], reverse=True)
    #print (words_freq[:100])

    #Generating wordcloud and saving as jpg image
    words_dict = dict(words_freq)
    WC_height = 1000
    WC_width = 1500
    WC_max_words = 200
    wordCloud = WordCloud(max_words=WC_max_words, height=WC_height, width=WC_width,stopwords=new_stopwords)
    wordCloud.generate_from_frequencies(words_dict)
    plt.title('Most frequently occurring bigrams connected by same colour and font size')
    plt.imshow(wordCloud, interpolation='bilinear')
    plt.axis("off")
    #plt.show()
    return wordCloud.to_file(filename_out)