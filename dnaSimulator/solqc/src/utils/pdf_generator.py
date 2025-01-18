from reportlab.pdfgen import canvas

from dnaSimulator.solqc.src.utils.content import Content

WIDTH = 595.25
HEIGHT = 841.89
TOP_MARGIN = 50
BOTTOM_MARGIN = 50
RIGHT_MARGIN = 30
TEXT_HEIGHT = 30
HEADLINE_HEIGHT = 2 * 30
IMAGE_WIDTH = 420
IMAGE_HEIGHT = IMAGE_WIDTH * 0.75

DEBUG = False

# TODO - Add functionality for regular text vs headlines
# TODO - Add functionality for text with new lines.
# TODO - Add functionality for quering images so we'll know how much to advance in y.

####   Utilities
def get_image_size(image_path):
    from PIL import Image
    im = Image.open(image_path)
    width, height = im.size
    return width, height

def get_image_ratio(image_path):
    width, height = get_image_size(image_path)
    return width / height


class PDFGenerator(object):
    def __init__(self, name):
        self.context = canvas.Canvas(name)
        self.name = name
        self.y = HEIGHT - TOP_MARGIN

        self.context.setFont("Courier", 14)

    def __enter__(self):
        if DEBUG:
            print("Object enter method has been called.")

        self.write_headline("Analysis Pipeline")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if DEBUG:
            print("Object exit method has been called.")
        self.close()

    def set_default_font(self):
        self.set_font('Courier', 14)

    def set_font(self, font_name, font_size):
        self.context.setFont(font_name, font_size)

    def write_headline(self, text):
        self.set_font('Courier', 20)
        self.context.drawCentredString(RIGHT_MARGIN + WIDTH/2, self.y, text)
        self.set_default_font()

        self.y -= HEADLINE_HEIGHT
        # TODO Implement

    def add_content_array(self, content_array, add_page_break=True):
        for content in content_array:
            self.add_content(content)

        if add_page_break:
            self.start_new_page()

    def add_content(self, content):
        content_type = content.type
        if content_type == Content.Type.TEXT:
            self.write_text(content())
        if content_type == Content.Type.IMAGE:
            self.add_image(content())
        if content_type == Content.Type.NEW_PAGE:
            self.start_new_page()

    def write_text(self, text):
        # Check if content should be in a new page.
        self.check_new_page(TEXT_HEIGHT)

        # Add content
        self.context.drawString(RIGHT_MARGIN, self.y, text)
        self.y -= TEXT_HEIGHT

    def add_image(self, image):
        # Check if content should be in a new page.
        self.check_new_page(IMAGE_HEIGHT)

        image_ratio = get_image_ratio(image)
        image_height = IMAGE_WIDTH / image_ratio
        # Add content
        self.context.drawImage(image, (WIDTH - IMAGE_WIDTH)/2, self.y - image_height, IMAGE_WIDTH, image_height)
        self.y -= IMAGE_HEIGHT

    def check_new_page(self, content_size):
        if self.y - content_size < BOTTOM_MARGIN:
            self.start_new_page()

    def start_new_page(self):
        self.context.showPage()
        self.set_default_font()
        self.y = HEIGHT - TOP_MARGIN

    def close(self):
        self.context.save()
        print (" - - {} file generated".format(self.name))
        if DEBUG:
            print("Pdf file closed")


if __name__ == "__main__":
    pdf = PDFGenerator("Sample.pdf")
    pdf.write_text("Testing pdf generator 2nd row")
    pdf.add_image('/Users/user/Dev/Academic/biology/SOLQC/temp/inputWIthXorLessError_withLong.png')
    pdf.add_image('/Users/user/Dev/Academic/biology/SOLQC/temp/reads_length_histogram.png')
    pdf.add_image('/Users/user/Dev/Academic/biology/SOLQC/temp/stacked_letters.png')
    pdf.add_image('/Users/user/Dev/Academic/biology/SOLQC/temp/stacked_letters_by_reads.png')
    pdf.add_image('/Users/user/Dev/Academic/biology/SOLQC/temp/SymbolDependent.png')
    pdf.add_image('/Users/user/Dev/Academic/biology/SOLQC/temp/variant_distribution.png')
    pdf.close()